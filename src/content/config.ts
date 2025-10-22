import { defineCollection, z } from 'astro:content';

const publications = defineCollection({
  type: 'content',
  schema: z.object({
    title: z.string(),
    authors: z.array(z.string()),
    journal: z.string(),
    year: z.number(),
    doi: z.string().optional(),
    pdf: z.string().optional(),
    preprint: z.string().optional(),
    abstract: z.string(),
    tags: z.array(z.string()),
    featured: z.boolean().default(false),
    citationCount: z.number().optional(),
    impactFactor: z.number().optional(),
    bibtex: z.string().optional(),
  }),
});

const projects = defineCollection({
  type: 'content',
  schema: z.object({
    title: z.string(),
    description: z.string(),
    status: z.enum(['ongoing', 'planned', 'completed']),
    startDate: z.date(),
    endDate: z.date().optional(),
    collaborators: z.array(z.string()).optional(),
    funding: z.string().optional(),
    tags: z.array(z.string()),
    image: z.string().optional(),
    repo: z.string().optional(),
    featured: z.boolean().default(false),
  }),
});

const blog = defineCollection({
  type: 'content',
  schema: z.object({
    title: z.string(),
    description: z.string(),
    pubDate: z.date(),
    category: z.enum(['python', 'swan', 'ww3', 'data-analysis', 'research-update', 'tutorial']),
    tags: z.array(z.string()),
    image: z.string().optional(),
    featured: z.boolean().default(false),
  }),
});

export const collections = {
  publications,
  projects,
  blog,
};
