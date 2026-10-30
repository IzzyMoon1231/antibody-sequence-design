# Reading design and evidence boundary

- Reader and updater routes are separate. A compact setup section selects a spreadsheet-style CSV; JSON is not a novice prerequisite.
- The opening overview answers availability and identity questions before interpretation. All supplied missing-artifact warnings remain visible.
- Each case groups humanness and structure. The offline mobile view collapses those groups, not the warnings. All plots remain available in the notebook and GitHub-friendly Markdown.
- Every plot has a question, source/cohort, axis guidance, computed median/count readout and limitation. Captions and mobile pages are generated from the same loaded run as the executed notebook. Scientific parsing, values, cohorts and figure geometry are unchanged.
- We use short, segmented narrative and expected-output explanations, rather than a fixed prose-to-code ratio. Relevant teaching guidance: [Berkeley notebook best practices](https://curriculum-guide.datahub.berkeley.edu/workflows/bestpractices/) and [EPFL interactive textbooks](https://www.epfl.ch/education/educational-initiatives/cede/teaching-interactively/jupyter-notebooks-for-education/teaching-and-learning-with-jupyter-notebooks/online-interactive-textbooks-with-jupyter-notebooks/).
- No renderer-specific hide-input promise is made. GitHub does not host raw HTML; `READ_RESULTS.md` uses relative image links and the notebook retains embedded PNGs. The mobile file is a local reading artifact, not a deployed app.

Rebuild the notebook source with `python make_notebook.py`, then restart and run it in the notebook folder. This replaces the notebook source, so preserve any saved execution you need first. The notebook's export step writes new reading views alongside new tables and figures. No prior export is overwritten.
