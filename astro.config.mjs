// @ts-check
import { defineConfig } from 'astro/config';
import mdx from '@astrojs/mdx';
import sitemap from '@astrojs/sitemap';
import tailwindcss from '@tailwindcss/vite';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import { readdirSync } from 'node:fs';

// Blog posts used to live at /blog/<file>.md/ ; keep those links working
const oldPostUrls = Object.fromEntries(
  readdirSync('./src/content/blog')
    .filter((f) => f.endsWith('.md'))
    .map((f) => [`/blog/${f}`, `/blog/${f.replace(/\.md$/, '')}/`])
);

// https://astro.build/config
export default defineConfig({
  site: 'https://dspelaez.github.io',
  i18n: {
    locales: ['en', 'es', 'fr'],
    defaultLocale: 'en',
    routing: { prefixDefaultLocale: false },
  },
  redirects: oldPostUrls,
  markdown: {
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeKatex],
  },
  integrations: [
    mdx({
      remarkPlugins: [remarkMath],
      rehypePlugins: [rehypeKatex],
    }),
    sitemap({
      i18n: { defaultLocale: 'en', locales: { en: 'en', es: 'es', fr: 'fr' } },
    }),
  ],
  vite: {
    plugins: [tailwindcss()]
  }
});