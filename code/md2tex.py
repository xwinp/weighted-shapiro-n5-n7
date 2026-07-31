#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert paper/n57_paper.md -> LaTeX for xelatex (xeCJK for Chinese)."""
import re, sys, html

SRC = "paper/n57_paper.md"
DST = "paper/n57_paper.tex"

with open(SRC, encoding="utf-8") as f:
    lines = f.readlines()

def esc_prose(s):
    """Escape LaTeX special chars in prose (math/code already extracted)."""
    s = s.replace("\\", r"\textbackslash{}")
    s = s.replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
    s = s.replace("_", r"\_").replace("{", r"\{").replace("}", r"\}")
    s = s.replace("~", r"\textasciitilde{}").replace("^", r"\textasciicircum{}")
    return s

def conv_inline(s):
    """Convert inline markdown in a prose span: extract math/code, escape, bold/italic."""
    # extract display $$...$$ (shouldn't appear inline but be safe) and inline $...$
    math = []
    def stash(m):
        math.append(m.group(0))
        return "\x00M%d\x00" % (len(math)-1)
    s = re.sub(r"\$\$.*?\$\$", stash, s, flags=re.S)
    s = re.sub(r"\$[^$]*\$", stash, s)
    # extract `code`
    code = []
    def stashc(m):
        code.append(m.group(1))
        return "\x00C%d\x00" % (len(code)-1)
    s = re.sub(r"`([^`]*)`", stashc, s)
    # escape prose
    s = esc_prose(s)
    # bold **...**  (may contain math placeholders)
    s = re.sub(r"\*\*(.+?)\*\*", lambda m: r"\textbf{"+m.group(1)+r"}", s, flags=re.S)
    # italic *...*
    s = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", lambda m: r"\textit{"+m.group(1)+r"}", s, flags=re.S)
    # restore code
    def restore_c(m):
        return r"\texttt{"+esc_prose(code[int(m.group(1))])+r"}"
    s = re.sub(r"\x00C(\d+)\x00", restore_c, s)
    # restore math (raw, unescaped)
    s = re.sub(r"\x00M(\d+)\x00", lambda m: math[int(m.group(1))], s)
    return s

out = []
in_itemize = False
in_enumerate = False
in_table = False
table_rows = []

def close_lists():
    global in_itemize, in_enumerate
    if in_itemize:
        out.append(r"\end{itemize}")
        in_itemize = False
    if in_enumerate:
        out.append(r"\end{enumerate}")
        in_enumerate = False

def flush_table():
    global in_table, table_rows
    if not in_table:
        return
    # table_rows: first row header, second row separator (---), rest body
    header = table_rows[0]
    body = table_rows[2:] if len(table_rows) > 2 else []
    ncol = len(header)
    spec = "l"*ncol
    out.append(r"\begin{tabular}{"+spec+r"}")
    out.append(r"\hline")
    out.append(" & ".join(conv_inline(c.strip()) for c in header) + r" \\")
    out.append(r"\hline")
    for row in body:
        out.append(" & ".join(conv_inline(c.strip()) for c in row) + r" \\")
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append("")
    in_table = False
    table_rows = []

# title from line 0
title_line = lines[0].rstrip("\n")[2:]
title_tex = conv_inline(title_line)

i = 1
# preamble lines until Abstract
body_start = None
for idx, ln in enumerate(lines[1:], start=1):
    if ln.startswith("## "):
        body_start = idx
        break
    # author / keywords / msc lines (author line is in \maketitle; skip it)
    s = ln.rstrip("\n")
    if s.strip() == "":
        continue
    if s.strip().startswith("**") and "@" in s:
        continue  # author line, shown by \maketitle
    out.append(conv_inline(s))
    out.append("")
out.append(r"\tableofcontents")
out.append(r"\bigskip")

i = body_start
while i < len(lines):
    raw = lines[i].rstrip("\n")
    raw = re.sub(r"^>\s?", "", raw)   # strip markdown blockquote marker (> Lemma/Proof)
    s = raw.strip()
    # display math block: single-line $$...$$ or multi-line $$ ... $$
    if s.startswith("$$"):
        close_lists(); flush_table()
        if s.endswith("$$") and len(s) > 4:
            body = s[2:-2]
            i += 1
        else:
            # collect until closing $$
            parts = [s[2:]]
            i += 1
            while i < len(lines):
                t = lines[i].rstrip("\n")
                if t.strip().endswith("$$"):
                    parts.append(t.strip()[:-2])
                    i += 1
                    break
                parts.append(t)
                i += 1
            body = "\n".join(parts)
        # long/wide equations: scale to text width (strip redundant size cmds)
        longq = len(body) > 160 or any(cmd in body for cmd in (r"\small", r"\scriptsize", r"\tiny", r"\footnotesize"))
        if longq:
            b2 = re.sub(r"\\(small|scriptsize|tiny|footnotesize)\s*", "", body)
            out.append(r"\begin{equation*}\resizebox{\textwidth}{!}{$\displaystyle "+b2+r"$}\end{equation*}")
        else:
            out.append("\\["+body+"\\]")
        out.append("")
        continue
    # table row
    if s.startswith("|") and s.endswith("|"):
        close_lists()
        cells = [c for c in s[1:-1].split("|")]
        if re.fullmatch(r"[\s:\-|]+", s):
            if not in_table:
                in_table = True
                table_rows = []
            table_rows.append(cells)
            i += 1
            continue
        if not in_table:
            in_table = True
            table_rows = []
        table_rows.append(cells)
        i += 1
        continue
    else:
        flush_table()
    # blank
    if s == "":
        close_lists()
        out.append("")
        i += 1
        continue
    # headers
    if s.startswith("### "):
        close_lists()
        t = s[4:]
        t = re.sub(r"^\d+\.\d+\s+", "", t)  # strip manual "3.1 " prefix
        out.append(r"\subsection{"+conv_inline(t)+"}")
        i += 1
        continue
    if s.startswith("## "):
        close_lists()
        t = s[3:]
        t = re.sub(r"^\d+\.\s+", "", t)  # strip manual "1. " prefix
        if t in ("Abstract", "References"):
            out.append(r"\addcontentsline{toc}{section}{"+conv_inline(t)+"}")
            out.append(r"\section*{"+conv_inline(t)+"}")
        else:
            out.append(r"\section{"+conv_inline(t)+"}")
        i += 1
        continue
    if s.startswith("# "):
        close_lists()
        out.append(r"\section*{"+conv_inline(s[2:])+"}")
        i += 1
        continue
    # bullet list
    if s.startswith("- "):
        if not in_itemize:
            close_lists()
            out.append(r"\begin{itemize}")
            in_itemize = True
        out.append(r"\item "+conv_inline(s[2:]))
        i += 1
        continue
    # numbered list  1. 2.
    if re.match(r"^\d+\.\s", s):
        if not in_enumerate:
            close_lists()
            out.append(r"\begin{enumerate}")
            in_enumerate = True
        out.append(r"\item "+conv_inline(re.sub(r"^\d+\.\s", "", s)))
        i += 1
        continue
    # horizontal rule ---
    if s == "---":
        close_lists()
        out.append(r"\par\noindent\rule{\textwidth}{0.4pt}")
        i += 1
        continue
    # plain paragraph
    close_lists()
    out.append(conv_inline(s))
    out.append("")
    i += 1

close_lists()
flush_table()

preamble = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=2.4cm]{geometry}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage{xeCJK}
\IfFontExistsTF{Noto Serif CJK SC}{\setCJKmainfont{Noto Serif CJK SC}}{\setCJKmainfont{SimSun}}
\IfFontExistsTF{Noto Sans CJK SC}{\setCJKsansfont{Noto Sans CJK SC}}{\setCJKsansfont{Microsoft YaHei}}
\usepackage{booktabs,array}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue!60!black,urlcolor=blue!60!black}
\newtheorem{theorem}{Theorem}
\newtheorem{proposition}{Proposition}
\newtheorem{lemma}{Lemma}
\newtheorem{corollary}{Corollary}
\theoremstyle{definition}
\newtheorem{definition}{Definition}
\theoremstyle{remark}
\newtheorem*{remark}{Remark}
\title{TITLE}
\author{AUTH}
\date{\today}
\begin{document}
\maketitle
\sloppy
\emergencystretch=3em
"""
# author line is lines[2] -> **薛炜鹏 (Weipeng Xue)** — 中山大学 ... ; conv_inline converts **..**
auth_tex = conv_inline(lines[2].rstrip("\n"))
preamble = preamble.replace("TITLE", title_tex.replace("{",r"\{").replace("}",r"\}"))
preamble = preamble.replace("AUTH", auth_tex)

tex = preamble + "\n".join(out) + "\n\\end{document}\n"
with open(DST, "w", encoding="utf-8") as f:
    f.write(tex)
print("wrote", DST, len(tex), "chars")
