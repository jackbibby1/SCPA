---
name: Error encountered
about: Any errors whilst running SCPA
title: ''
labels: ''
assignees: ''

---

**The specific error you're running into**
Either copy/paste the error message, or provide a screenshot

**To Reproduce**
The code to reproduce the error e.g.

```
library(SCPA); library(Seurat); library(msigdbr)

p1 <- seurat_extract(df, meta1 = "hour", value_meta1 = "0")
p2 <- seurat_extract(df, meta1 = "hour", value_meta1 = "24")

pathways <- msigdbr("Homo sapiens", "H") %>% format_pathways()

scpa_out <- compare_pathways(list(p1, p2), pathways)
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Additional context**
Add any other context about the problem here.
