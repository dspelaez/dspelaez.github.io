# Academic CV Website - Daniel Peláez-Zapata

A modern academic website built with Astro, showcasing research, publications, projects, and tutorials in physical oceanography.

## Features

- 🎨 Modern, clean design with dark mode support
- 📱 Fully responsive and mobile-friendly
- 📝 Markdown-based content management
- 🔍 SEO optimized with Schema.org markup
- 📊 Publication tracking with citations and impact factors
- 🚀 Static site generation (fast loading)
- 📖 Math rendering support with KaTeX
- 💻 Code syntax highlighting
- 🔗 Social links integration (GitHub, Google Scholar, ORCID)

## Quick Start

### Prerequisites

- Node.js 18+ and npm

### Installation

```bash
# Install dependencies
npm install

# Start development server
npm run dev

# Build for production
npm run build

# Preview production build
npm run preview
```

## Content Management

All content is managed through Markdown files in the `src/content` directory:

### Publications

Add new publications to `src/content/publications/`:

```markdown
---
title: "Your Paper Title"
authors: ["Author 1", "Author 2"]
journal: "Journal Name"
year: 2024
doi: "10.xxxx/xxxxx"
pdf: "/papers/yourpaper.pdf"
abstract: "Your abstract here..."
tags: ["tag1", "tag2"]
citationCount: 10
impactFactor: 3.5
bibtex: |
  @article{...}
---
```

### Projects

Add projects to `src/content/projects/`:

```markdown
---
title: "Project Title"
description: "Brief description"
status: "ongoing" # or "planned" or "completed"
startDate: 2024-01-01
tags: ["tag1", "tag2"]
repo: "https://github.com/..."
---

## Details

Project details in markdown...
```

### Blog Posts

Add blog posts/tutorials to `src/content/blog/`:

```markdown
---
title: "Post Title"
description: "Brief description"
pubDate: 2024-01-15
category: "python" # python, swan, ww3, data-analysis, research-update, tutorial
tags: ["tag1", "tag2"]
---

Your content here with code blocks, math, etc...
```

## Customization

### Personal Information

Update the following files:

1. **About Section**: Edit `src/components/sections/About.astro`
   - Update bio, education, and skills

2. **Contact Information**: Edit `src/components/Footer.astro` and `src/components/sections/Contact.astro`
   - Update email, GitHub, Google Scholar, ORCID links

3. **Profile Photo**: Replace `public/profile.jpg` with your photo

4. **Site Metadata**: Edit `astro.config.mjs`
   - Update `site` URL

### Styling

- Global styles: `src/styles/global.css`
- Tailwind configuration: Uses Tailwind CSS 4.x

## Deployment to GitHub Pages

### Setup

1. Push your code to a GitHub repository
2. Go to repository Settings > Pages
3. Set Source to "GitHub Actions"
4. The site will automatically deploy on every push to `main`

### Custom Domain (Optional)

1. Add a `CNAME` file to the `public` directory with your domain
2. Configure DNS settings with your domain provider

## Project Structure

```
/
├── public/              # Static assets
│   ├── profile.jpg      # Your profile photo
│   └── favicon.png      # Site icon
├── src/
│   ├── components/      # Reusable components
│   │   ├── sections/    # Homepage sections
│   │   ├── Header.astro
│   │   └── Footer.astro
│   ├── content/         # Markdown content
│   │   ├── publications/
│   │   ├── projects/
│   │   └── blog/
│   ├── layouts/         # Page layouts
│   ├── pages/           # Site pages
│   │   ├── index.astro
│   │   ├── publications/
│   │   ├── projects/
│   │   └── blog/
│   └── styles/          # Global styles
├── astro.config.mjs     # Astro configuration
└── package.json
```

## Tech Stack

- **Framework**: Astro 4.x
- **Styling**: Tailwind CSS 4.x
- **Content**: Markdown + MDX
- **Math**: KaTeX
- **Icons**: Heroicons (via inline SVG)
- **Deployment**: GitHub Pages

## Development Commands

| Command                | Action                                           |
|:-----------------------|:-------------------------------------------------|
| `npm install`          | Install dependencies                             |
| `npm run dev`          | Start dev server at `localhost:4321`             |
| `npm run build`        | Build production site to `./dist/`               |
| `npm run preview`      | Preview production build locally                 |

## Contributing

This is a personal academic website. Feel free to fork and adapt for your own use.

## License

MIT License - feel free to use this template for your own academic website.

## Support

For issues or questions:
- Check the [Astro documentation](https://docs.astro.build)
- Open an issue on GitHub

---

Built with [Astro](https://astro.build) and [Tailwind CSS](https://tailwindcss.com)
