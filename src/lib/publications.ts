import { readFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import type { Publication } from '../components/PublicationCard.astro';

let cache: Publication[] | null = null;

/** My publications (the ORCID/CrossRef cache written by the Publications section). */
export async function getMyPublications(): Promise<Publication[]> {
  if (!cache) {
    try {
      const data = JSON.parse(await readFile(resolve('./public/crossref-publications.json'), 'utf-8'));
      cache = data.publications || [];
    } catch {
      cache = [];
    }
  }
  return cache!;
}

/** Find my publications by DOI (or by URL for works without DOI), in the given order. */
export async function findMyPapers(ids: string[]): Promise<Publication[]> {
  const pubs = await getMyPublications();
  const norm = (s: string) => s.toLowerCase().replace(/^https?:\/\/(dx\.)?doi\.org\//, '').trim();
  return ids
    .map((id) => pubs.find((p) => (p.doi && norm(p.doi) === norm(id)) || (p.url && p.url === id)))
    .filter((p): p is Publication => Boolean(p));
}
