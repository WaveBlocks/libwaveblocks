SVG2PDF=inkscape -z -T -A

.PHONY all: manual.pdf

.PHONY clean:
	rm -f *.pdf
	rm -f *.log
	rm -f *.pdf_tex
	rm -f *.aux
	rm -f *.pdf~
	rm -f *.tex~
	rm -f *.svg~
	rm -f *.toc

# DOT files (you need graphviz to compile these files)

hawp_inheritance.svg~: hawp_inheritance.dot
	dot -Tsvg $< -o $@

hawp_inheritance.pdf: hawp_inheritance.svg~
	$(SVG2PDF) $@ $<

# SVG Files

shape_slicing.pdf: shape_slicing.svg
	$(SVG2PDF) $@ $<

shape_extension.pdf: shape_extension.svg
	$(SVG2PDF) $@ $<

shape_example.pdf: shape_example.svg
	$(SVG2PDF) $@ $<

shape_enumerator.pdf: shape_enumerator.svg
	$(SVG2PDF) $@ $<

grad_scatter_stencil.pdf: grad_scatter_stencil.svg
	$(SVG2PDF) $@ $<

grad_gather_stencil.pdf: grad_gather_stencil.svg
	$(SVG2PDF) $@ $<

basis_eval_stencil.pdf: basis_eval_stencil.svg
	$(SVG2PDF) $@ $<

manual.pdf: manual.tex header_math.tex header_code.tex basis_eval_stencil.pdf grad_gather_stencil.pdf grad_scatter_stencil.pdf shape_enumerator.pdf shape_example.pdf shape_extension.pdf shape_slicing.pdf hawp_inheritance.pdf
	pdflatex manual.tex
#	pdflatex --shell-escape manual.tex	

