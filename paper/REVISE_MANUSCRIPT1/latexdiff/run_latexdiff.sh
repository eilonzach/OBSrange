#!/bin/bash

draft="OBSrange.tex"
revision="../OBSrange.tex"

latexdiff ${draft} ${revision} --flatten > diff.tex

# latexdiffcite file ${draft} ${revision}

