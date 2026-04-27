## About

This repository contains the source code for my personal website, 
built with [Hugo](https://gohugo.io/) — a fast, open-source static site generator.
The site serves as a space for personal essays, notes, and meanderings across topics I find genuinely interesting.

---

## Tech Stack

| Layer | Tool |
|---|---|
| Static Site Generator | [Hugo](https://gohugo.io/) |
| Markup | Markdown (via Goldmark) |
| Math Rendering | LaTeX (via Goldmark passthrough — supports `$$...$$`, `\[...\]`, `\(...\)`) |
| Styling | Custom CSS |
| Deployment | GitHub Actions → [karso.xyz](https://karso.xyz) |

---

## Project Structure

```
karsosp/
├── .github/workflows/   # CI/CD pipeline (auto-deploy on push)
├── archetypes/          # Hugo content templates
├── assets/css/          # Custom stylesheets
├── content/             # Site content (Markdown posts/pages)
├── layouts/             # HTML templates overriding theme defaults
├── static/              # Static assets (images, fonts, etc.)
└── hugo.toml            # Site configuration
```

---

This repository is personal and not intended for redistribution. Feel free to browse the source for reference.
