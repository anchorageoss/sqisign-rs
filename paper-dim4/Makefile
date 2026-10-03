# Regenerates every figure and table of the paper from the repository's
# documents (BENCH.md, docs/COMPACT_R3.md, COMPRESSION.md, README.md).
# `make figures` writes figures/data/*.csv and *.tex (read by pgfplots and by
# the macros in main.tex), figures/out/*.pdf|svg|png (matplotlib renderings
# of the data figure) and tables/*.tex; then figures/bundle.py inlines all
# of it into main_standalone.tex, the one file to paste into Overleaf next
# to refs.bib.  `make check` runs the generation twice and verifies the
# output is byte-identical.  Same layout as prism-rs paper/prism-suf.
PY ?= python3
FIGS := values fig02_tradeoff
TABS := tab1_parameters tab2_sizes tab3_bench tab4_ms tab5_loop

.PHONY: figures check clean pdf words
figures:
	@for s in $(FIGS) $(TABS); do SOURCE_DATE_EPOCH=0 $(PY) figures/$$s.py || exit 1; done
	@$(PY) figures/bundle.py

check: figures
	@cd figures && find data out ../tables ../main_standalone.tex -type f | sort | xargs sha256sum > /tmp/dim4_paper_figs_1.txt
	@$(MAKE) --no-print-directory figures
	@cd figures && find data out ../tables ../main_standalone.tex -type f | sort | xargs sha256sum > /tmp/dim4_paper_figs_2.txt
	@cmp /tmp/dim4_paper_figs_1.txt /tmp/dim4_paper_figs_2.txt && echo "byte-identical on the second run"

clean:
	rm -f figures/data/* figures/out/* tables/*.tex main_standalone.tex

# Local render (TeX Live + bibtex, three passes); the author's canonical
# build is Overleaf.  Output: main_standalone.pdf, intermediates in build/.
pdf: figures
	@mkdir -p build && cp main_standalone.tex refs.bib build/
	@cd build && pdflatex -interaction=nonstopmode main_standalone.tex > p1.log 2>&1 || (grep -n '^!' -A3 p1.log; exit 1)
	@cd build && bibtex main_standalone > b.log 2>&1 || true
	@cd build && pdflatex -interaction=nonstopmode main_standalone.tex > p2.log 2>&1 && pdflatex -interaction=nonstopmode main_standalone.tex > p3.log 2>&1
	@cp build/main_standalone.pdf main_standalone.pdf
	@echo "errors: $$(grep -c '^!' build/p3.log)"; grep -n 'Overfull\|pgfplots Warning\|Citation.*undefined\|Reference.*undefined' build/main_standalone.log | head; echo "pages: $$(pdfinfo main_standalone.pdf | awk '/^Pages/{print $$2}')"

# The ticket's register checks: banned words and em dashes.
words:
	@! grep -v "^[[:space:]]*%" main.tex | grep -n -i -E "\bnovel\b|\bsignificant|\brobust|\belegant|we believe|\bclearly\b" || (echo "banned word above"; exit 1)
	@! grep -v "^[[:space:]]*%" main.tex | LC_ALL=C.UTF-8 grep -n -P "\x{2014}|---" || (echo "em dash above"; exit 1)
	@echo "register: clean"
