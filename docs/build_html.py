#!/usr/bin/env python3
"""Render the research notes (and the plan documents they point to) as self-contained HTML pages.

    uv run --no-project --with markdown-it-py==3.0.0 python docs/build_html.py

Writes docs/html/. Every page carries its own CSS, no scripts and no fonts from the network, so a single
file can be opened, emailed or printed. The markdown stays the source: rebuild after editing it.
"""

from __future__ import annotations

import html
import re
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from markdown_it import MarkdownIt

DOCS = Path(__file__).resolve().parent
OUT = DOCS / "html"


@dataclass(frozen=True)
class Doc:
    src: str  # relative to docs/
    label: str  # the number the file carries, or a role
    title: str
    blurb: str
    group: str  # "research" or "plan"

    @property
    def out(self) -> str:
        return Path(self.src).stem.lower() + ".html"


ALL = [
    Doc("research/SUMMARY.md", "Start", "Summary", "Top recommendations, likely target tiers, what is still open and what to ask the organisers.", "research"),
    Doc("research/01_binder_filtering.md", "01", "Filtering binders", "Which in-silico metrics predicted wet-lab binding in Adaptyv data, with recommended thresholds.", "research"),
    Doc("research/02_target_candidates.md", "02", "Target candidates", "Breast-cancer-relevant extracellular proteins (HER2, EGFR, HER3, TROP2, Nectin-4, PD-L1, FGFR2) in likelihood tiers.", "research"),
    Doc("research/03_generation_tools.md", "03", "Generation tools", "Binder design tools compared, with a recommended setup for one 10-hour day.", "research"),
    Doc("research/04_adaptyv_rules_and_assay.md", "04", "Adaptyv rules and assay", "Submission rules per round, how the assay and a hit are defined, and a checklist of sequence constraints.", "research"),
    Doc("research/05_winners_and_jury.md", "05", "Winners and jury", "Past winners and jury criteria across events, ending in a pitch checklist.", "research"),
    Doc("research/06_track_a_antibody_optimisation.md", "06", "Track A antibody optimisation", "Improving an existing breast cancer scFv in one blind round: methods, positions, known mutations, ranking and a one-day plan.", "research"),
    Doc("research/07_optimised_binder_data.md", "07", "Optimised binder data", "What the public Adaptyv data show about optimised and mutated binders compared with their parents.", "research"),
    Doc("filtering_strategy.md", "Plan", "Filtering strategy", "The single filter specification: stages, thresholds, and what changes when the metric is announced.", "plan"),
    Doc("challenge_notes.md", "Plan", "Challenge notes", "What the event page says, precedents for the format, the likely challenge and questions for the organisers.", "plan"),
    Doc("track_a_glossary.md", "Plan", "Track A glossary", "Plain definitions of the antibody, assay and scoring terms used in the Track A notes.", "plan"),
]
BY_NAME = {Path(d.src).name: d for d in ALL}

CSS = """
:root{
  --bg:#f5f7fa; --surface:#ffffff; --fg:#19212e; --muted:#576274; --line:#d9dfe8;
  --accent:#1c5d89; --accent-soft:#e3eef7; --row:#eff3f8; --code-bg:#e9edf3;
  --flag-bg:#fde9c4; --flag-fg:#6b4000; --chip-on:#ffffff;
  --font:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;
  --mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --bg:#0e131a; --surface:#141a23; --fg:#e3e8ef; --muted:#9ba6b6; --line:#28313e;
  --accent:#74b4e3; --accent-soft:#15293a; --row:#1a212c; --code-bg:#1d2531;
  --flag-bg:#3b2a08; --flag-fg:#f1c675; --chip-on:#0e131a; color-scheme:dark}}
:root[data-theme="dark"]{
  --bg:#0e131a; --surface:#141a23; --fg:#e3e8ef; --muted:#9ba6b6; --line:#28313e;
  --accent:#74b4e3; --accent-soft:#15293a; --row:#1a212c; --code-bg:#1d2531;
  --flag-bg:#3b2a08; --flag-fg:#f1c675; --chip-on:#0e131a; color-scheme:dark}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
@media (prefers-reduced-motion:no-preference){html{scroll-behavior:smooth}}
body{margin:0;background:var(--bg);color:var(--fg);font:17px/1.65 var(--font);overflow-wrap:anywhere;
  padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px));padding-block:0 4rem}
a{color:var(--accent);text-underline-offset:.15em}
a:focus-visible,summary:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:3px}
.wrap{max-width:78rem;margin:0 auto}
.topbar{display:flex;flex-wrap:wrap;gap:.25rem 1.5rem;justify-content:space-between;align-items:baseline;padding-block:1.1rem;font-size:.9rem;color:var(--muted)}
.topbar a{color:var(--muted);text-decoration:none}.topbar a:hover{color:var(--accent);text-decoration:underline}
.topbar strong{color:var(--fg);font-weight:650;letter-spacing:.01em}
.layout{display:grid;grid-template-columns:15.5rem minmax(0,1fr);gap:3.5rem;align-items:start}
.toc{font-size:.88rem;line-height:1.4}
.toc-side{position:sticky;top:1.25rem;max-height:calc(100vh - 2.5rem);overflow:auto;padding-right:.5rem}
.toc-title{margin:0 0 .6rem;font-size:.74rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);font-weight:650}
.toc ol{list-style:none;margin:0;padding:0}
.toc li{margin:0}.toc li a{display:block;padding:.28rem 0;color:var(--muted);text-decoration:none}
.toc li a:hover{color:var(--accent)}
.toc li.sub a{padding-left:.9rem;font-size:.84rem}
.toc-top{display:none;margin:0 0 1.5rem;border:1px solid var(--line);border-radius:8px;background:var(--surface);padding:.2rem 1rem}
.toc-top summary{cursor:pointer;padding:.6rem 0;font-weight:600;font-size:.92rem}
@media (max-width:1000px){.layout{grid-template-columns:minmax(0,1fr);gap:0}.toc-side{display:none}.toc-top{display:block}}
.eyebrow{margin:0 0 .5rem;font-size:.78rem;letter-spacing:.09em;text-transform:uppercase;color:var(--accent);font-weight:650}
h1{margin:0 0 .6rem;font-size:clamp(1.7rem,3.2vw,2.35rem);line-height:1.18;letter-spacing:-.015em;text-wrap:balance;max-width:34ch}
.meta{display:flex;flex-wrap:wrap;gap:.2rem 1.25rem;margin:0 0 2rem;color:var(--muted);font-size:.88rem}
main{min-width:0}
main>:is(p,ul,ol,blockquote,h2,h3,h4,hr){max-width:70ch}
h2{margin:2.8rem 0 .8rem;font-size:1.45rem;line-height:1.25;letter-spacing:-.01em;text-wrap:balance;scroll-margin-top:1rem}
h3{margin:2rem 0 .5rem;font-size:1.12rem;line-height:1.3;scroll-margin-top:1rem}
h4{margin:1.5rem 0 .4rem;font-size:1rem}
.anchor{margin-left:.45rem;color:var(--muted);text-decoration:none;opacity:0;font-weight:400}
h2:hover .anchor,h3:hover .anchor,.anchor:focus{opacity:1}
p{margin:0 0 1rem}
ul,ol{margin:0 0 1.1rem;padding-left:1.35rem}li{margin:.4rem 0}li::marker{color:var(--accent)}
li>ul,li>ol{margin:.3rem 0 .2rem}
strong{font-weight:650}
hr{border:0;border-top:1px solid var(--line);margin:2.5rem 0}
blockquote{margin:0 0 1.1rem;padding:.1rem 0 .1rem 1.1rem;border-left:3px solid var(--accent);color:var(--muted)}
code{font:.86em/1.4 var(--mono);background:var(--code-bg);padding:.12em .38em;border-radius:4px;overflow-wrap:anywhere}
a code{color:inherit}
pre{overflow-x:auto;background:var(--code-bg);padding:1rem;border-radius:8px;font-size:.88rem}pre code{background:none;padding:0}
.table-wrap{overflow-x:auto;margin:1.3rem 0 1.6rem;background:var(--surface);border:1px solid var(--line);border-radius:8px}
table{border-collapse:collapse;width:100%;min-width:34rem;font-size:.9rem;line-height:1.45;font-variant-numeric:tabular-nums}
th{background:var(--row);text-align:left;font-weight:650;white-space:nowrap}
th,td{padding:.6rem .8rem;border-bottom:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}
td{min-width:5.5rem}
tr:last-child td{border-bottom:0}
tbody tr:hover td{background:var(--row)}
td:first-child{font-weight:560}
.chip{display:inline-block;padding:.05rem .6rem;border-radius:99px;font-size:.8rem;font-weight:650;white-space:nowrap;border:1px solid var(--accent);color:var(--accent)}
.chip.likely{background:var(--accent);color:var(--chip-on)}
.chip.long{border-color:var(--line);color:var(--muted)}
mark.flag{background:var(--flag-bg);color:var(--flag-fg);padding:.04em .4em;border-radius:4px;font-size:.78em;font-weight:700;letter-spacing:.03em;white-space:nowrap}
.cite{padding:0 .25em;border-radius:3px;background:var(--accent-soft);color:var(--accent);text-decoration:none;font-size:.82em;white-space:nowrap}
.cite:hover{text-decoration:underline}
.src-id{scroll-margin-top:1rem;font-weight:650}
.src-id:target{background:var(--flag-bg);color:var(--flag-fg);border-radius:3px;padding:0 .25em}
td.src-id:target{background:var(--flag-bg)}
.pager{display:flex;flex-wrap:wrap;justify-content:space-between;gap:1rem;margin:3.5rem 0 0;padding-top:1.5rem;border-top:1px solid var(--line);max-width:70ch;font-size:.95rem}
.pager small{display:block;color:var(--muted);font-size:.78rem;letter-spacing:.06em;text-transform:uppercase}
.pager a{text-decoration:none}.pager a:hover strong{text-decoration:underline}
.pager .next{margin-left:auto;text-align:right}
footer{margin:2.5rem 0 0;color:var(--muted);font-size:.82rem;max-width:70ch}
.intro{max-width:62ch;color:var(--muted);margin:0 0 2.2rem;font-size:1.05rem}
.group-title{margin:2.4rem 0 .6rem;font-size:.78rem;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);font-weight:650}
.doc-list{display:flex;flex-direction:column;border:1px solid var(--line);border-radius:10px;background:var(--surface);overflow:hidden}
.doc-row{display:grid;grid-template-columns:4.2rem minmax(0,1fr) auto;gap:1rem;align-items:baseline;padding:1rem 1.2rem;color:inherit;text-decoration:none;border-top:1px solid var(--line)}
.doc-row:first-child{border-top:0}.doc-row:hover{background:var(--row)}
.doc-row .no{font-size:.82rem;color:var(--accent);font-weight:700;letter-spacing:.04em;text-transform:uppercase}
.doc-row .t{display:block;font-weight:650;font-size:1.05rem;line-height:1.35}.doc-row:hover .t{color:var(--accent)}
.doc-row .b{display:block;color:var(--muted);font-size:.92rem;line-height:1.5;margin-top:.15rem}
.doc-row .m{font-size:.82rem;color:var(--muted);text-align:right;white-space:nowrap;display:flex;flex-direction:column;gap:.25rem;align-items:flex-end}
@media (max-width:640px){.doc-row{grid-template-columns:minmax(0,1fr);gap:.3rem}.doc-row .m{flex-direction:row;gap:.8rem;text-align:left}}
@media print{body{background:#fff;color:#000;font-size:11pt;padding:0}.topbar,.toc,.pager,.anchor{display:none}.layout{display:block}
  .table-wrap{overflow:visible;border-color:#999}table{min-width:0;font-size:9pt}a{color:inherit}}
"""

md = MarkdownIt("commonmark", {"breaks": True, "html": False}).enable("table").enable("strikethrough")
TAG = re.compile(r"(<[^>]+>)")


def slugify(text: str, used: set[str]) -> str:
    base = re.sub(r"[^a-z0-9]+", "-", re.sub(r"[`*_]", "", text).lower()).strip("-") or "section"
    slug, n = base, 2
    while slug in used:
        slug, n = f"{base}-{n}", n + 1
    used.add(slug)
    return slug


def map_text(markup: str, fn) -> str:
    """Apply fn to the text between tags, leaving links, code and marks alone."""
    out, skip = [], 0
    for part in TAG.split(markup):
        if part.startswith("<"):
            found = re.match(r"<(/?)(a|code|pre|mark)\b", part)
            if found:
                skip += -1 if found.group(1) else 1
            out.append(part)
        else:
            out.append(part if skip > 0 else fn(part))
    return "".join(out)


def decorate(body: str) -> str:
    """Anchors for source entries, links for citations, URLs and document paths, flags, chips, table wrappers."""
    seen: set[str] = set()

    def source_id(m: re.Match) -> str:
        if m.group(2) in seen:
            return m.group(0)
        seen.add(m.group(2))
        return f'{m.group(1)}<span id="{m.group(2)}" class="src-id">{m.group(2)}</span>'

    def cell_id(m: re.Match) -> str:
        if m.group(1) in seen:
            return m.group(0)
        seen.add(m.group(1))
        return f'<td id="{m.group(1)}" class="src-id">{m.group(1)}</td>'

    body = re.sub(r"(<p>|<br />\n?)(S\d+)(?= )", source_id, body)
    body = re.sub(r"<td>(S\d+)</td>", cell_id, body)

    # [S12], [S12, S4], [S24 to S27], [S31, my arithmetic]: each S-id becomes a link, the rest stays text
    cite = re.compile(r"\[(S\d+(?:(?:,\s*|\s+to\s+)S\d+)*)((?:,\s[^\]]*)?)\]")
    link_ids = lambda ids: re.sub(r"S\d+", lambda s: f'<a class="cite" href="#{s.group(0)}">{s.group(0)}</a>' if s.group(0) in seen else s.group(0), ids)
    body = map_text(body, lambda t: cite.sub(lambda m: "[" + link_ids(m.group(1)) + m.group(2) + "]", t))
    body = map_text(body, lambda t: re.sub(r"\bUNVERIFIED\b", '<mark class="flag">UNVERIFIED</mark>', t))

    def url(m: re.Match) -> str:
        link, trail = m.group(0), ""
        while link and link[-1] in ".,;:)]":
            link, trail = link[:-1], link[-1] + trail
        return f'<a href="{link}">{link}</a>{trail}'

    body = map_text(body, lambda t: re.sub(r'https?://[^\s<>"]+', url, t))
    body = re.sub(r'<a href="(https?://[^"]+)"', r'<a href="\1" target="_blank" rel="noopener"', body)

    def doc_path(m: re.Match) -> str:
        target = BY_NAME.get(Path(m.group(1)).name)
        return f'<a href="{target.out}"><code>{m.group(1)}</code></a>' if target else m.group(0)

    body = re.sub(r"<code>([\w./-]+\.md)</code>", doc_path, body)
    body = re.sub(r"<td>(Likely)</td>", r'<td><span class="chip likely">\1</span></td>', body)
    body = re.sub(r"<td>(Possible)</td>", r'<td><span class="chip">\1</span></td>', body)
    body = re.sub(r"<td>(Long tail)</td>", r'<td><span class="chip long">\1</span></td>', body)
    body = body.replace("<table>", '<div class="table-wrap"><table>').replace("</table>", "</table></div>")
    return re.sub(
        r'<h([23]) id="([^"]+)">(.*?)</h\1>',
        lambda m: f'<h{m[1]} id="{m[2]}">{m[3]}<a class="anchor" href="#{m[2]}" aria-label="Link to this section">#</a></h{m[1]}>',
        body,
    )


def page(title: str, body: str, *, crumb: str = "") -> str:
    return (
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        f"<title>{html.escape(title)}</title>\n<style>{CSS}</style>\n</head>\n<body>\n"
        f'<div class="wrap">\n<div class="topbar"><strong><a href="index.html">PROTECher research notes</a></strong><span>{crumb}</span></div>\n'
        f"{body}\n</div>\n</body>\n</html>\n"
    )


def build_doc(doc: Doc) -> tuple[str, dict]:
    text = (DOCS / doc.src).read_text()
    tokens = md.parse(text)
    used: set[str] = set()
    toc: list[tuple[int, str, str]] = []
    h1_html, drop = doc.title, None
    for i, tok in enumerate(tokens):
        if tok.type != "heading_open":
            continue
        inline = tokens[i + 1]
        level = int(tok.tag[1])
        if level == 1 and drop is None:
            h1_html, drop = md.renderInline(inline.content), i
            continue
        tok.attrSet("id", slugify(inline.content, used))
        if level in (2, 3):
            toc.append((level, re.sub(r"[`*]", "", inline.content), tok.attrGet("id")))
    if drop is not None:
        del tokens[drop : drop + 3]
    body = decorate(md.renderer.render(tokens, md.options, {}))

    words = len(re.findall(r"\w+", text))
    flags = text.count("UNVERIFIED")
    items = "".join(f'<li class="{"sub" if lvl == 3 else ""}"><a href="#{slug}">{html.escape(label)}</a></li>' for lvl, label, slug in toc)
    side = f'<nav class="toc toc-side" aria-label="On this page"><p class="toc-title">On this page</p><ol>{items}</ol></nav>'
    top = f'<details class="toc toc-top"><summary>On this page</summary><ol>{items}</ol></details>'

    same = [d for d in ALL if d.group == doc.group]
    at = same.index(doc)
    prev = f'<a class="prev" href="{same[at - 1].out}"><small>Previous</small><strong>{html.escape(same[at - 1].title)}</strong></a>' if at else ""
    nxt = f'<a class="next" href="{same[at + 1].out}"><small>Next</small><strong>{html.escape(same[at + 1].title)}</strong></a>' if at + 1 < len(same) else ""
    eyebrow = {"Start": "Research summary", "Plan": "Plan document"}.get(doc.label, f"Research note {doc.label}")
    meta = f"<span>{max(1, round(words / 220))} min read</span><span>{words:,} words</span>"
    if flags:
        meta += f'<span><mark class="flag">UNVERIFIED</mark> {flags} marked claim{"s" if flags != 1 else ""}</span>'
    footer = f"<footer>Rendered from <code>docs/{doc.src}</code> on {date.today().isoformat()}. The markdown file is the source; rebuild with <code>docs/build_html.py</code>.</footer>"
    content = (
        f'<div class="layout">{side}<main>{top}<p class="eyebrow">{html.escape(eyebrow)}</p><h1>{h1_html}</h1>'
        f'<p class="meta">{meta}</p>{body}<nav class="pager" aria-label="Neighbouring documents">{prev}{nxt}</nav>{footer}</main></div>'
    )
    plain_title = re.sub(r"<[^>]+>", "", h1_html)
    return page(f"{plain_title} | PROTECher research", content, crumb=f'<a href="index.html">All notes</a>'), {"words": words, "flags": flags}


def build_index(stats: dict[str, dict]) -> str:
    def rows(group: str) -> str:
        out = []
        for d in (d for d in ALL if d.group == group):
            s = stats[d.src]
            flag = f'<span><mark class="flag">UNVERIFIED</mark> {s["flags"]}</span>' if s["flags"] else ""
            out.append(
                f'<a class="doc-row" href="{d.out}"><span class="no">{d.label}</span>'
                f'<span><span class="t">{html.escape(d.title)}</span><span class="b">{html.escape(d.blurb)}</span></span>'
                f'<span class="m"><span>{max(1, round(s["words"] / 220))} min</span>{flag}</span></a>'
            )
        return '<div class="doc-list">' + "".join(out) + "</div>"

    body = (
        '<main style="padding-top:1rem"><p class="eyebrow">Breast Cancer Awareness Hackathon, Zürich, 31 October 2026</p>'
        "<h1>PROTECher research notes</h1>"
        '<p class="intro">Background for the Adaptyv protein design challenge, written 24 to 30 October. Read the summary first. '
        "A claim marked UNVERIFIED was not confirmed against a primary source, and the count beside each note says how many there are.</p>"
        f'<p class="group-title">Research</p>{rows("research")}<p class="group-title">Plan and glossary</p>{rows("plan")}'
        f"<footer>Built {date.today().isoformat()} from the markdown files in <code>docs/</code> with <code>docs/build_html.py</code>.</footer></main>"
    )
    return page("PROTECher research notes", body)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    stats = {}
    for doc in ALL:
        markup, stats[doc.src] = build_doc(doc)
        (OUT / doc.out).write_text(markup)
    (OUT / "index.html").write_text(build_index(stats))
    print(f"wrote {len(ALL) + 1} pages to {OUT}")


if __name__ == "__main__":
    main()
