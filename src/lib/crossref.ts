interface CrossrefPublication {
  title: string;
  authors: string;
  journal: string;
  year: number;
  doi: string;
  abstract: string;
  citationKey: string;
  itemType: string;
  citationCount?: number;
  volume?: string;
  issue?: string;
  pages?: string;
}

interface CrossrefAuthor {
  given?: string;
  family?: string;
  sequence?: string;
}

interface CrossrefWork {
  title?: string[];
  author?: CrossrefAuthor[];
  'container-title'?: string[];
  'short-container-title'?: string[];
  publisher?: string;
  published?: {
    'date-parts'?: number[][];
  };
  issued?: {
    'date-parts'?: number[][];
  };
  DOI?: string;
  abstract?: string;
  type?: string;
  'is-referenced-by-count'?: number;
  volume?: string;
  issue?: string;
  page?: string;
}

/**
 * Fetch publications from CrossRef API using ORCID
 * @param orcid - The ORCID identifier
 * @param email - Email for CrossRef polite pool (optional but recommended)
 * @returns Array of publications
 */
export async function fetchPublicationsFromCrossref(
  orcid: string,
  email: string = 'noreply@example.com'
): Promise<CrossrefPublication[]> {
  const url = `https://api.crossref.org/works?filter=orcid:${orcid}&rows=100&mailto=${email}&sort=published&order=desc`;

  try {
    console.log(`Fetching publications from CrossRef for ORCID: ${orcid}`);
    const response = await fetch(url);

    if (!response.ok) {
      throw new Error(`CrossRef API error: ${response.status} ${response.statusText}`);
    }

    const data = await response.json();
    const items = data.message?.items || [];

    console.log(`Found ${items.length} publications from CrossRef`);

    return items.map(transformWork).filter(Boolean) as CrossrefPublication[];
  } catch (error) {
    console.error('Failed to fetch from CrossRef:', error);
    throw error;
  }
}

/**
 * Transform CrossRef work data to our publication format
 */
function transformWork(work: CrossrefWork): CrossrefPublication | null {
  try {
    // Extract authors
    const authors = work.author
      ?.map((a: CrossrefAuthor) => {
        const given = a.given || '';
        const family = a.family || '';
        return `${given} ${family}`.trim();
      })
      .filter(name => name.length > 0)
      .join(', ') || 'Unknown';

    // Extract year from published or issued date
    const year =
      work.published?.['date-parts']?.[0]?.[0] ||
      work.issued?.['date-parts']?.[0]?.[0] ||
      new Date().getFullYear();

    // Extract journal name
    const journal =
      work['container-title']?.[0] ||
      work['short-container-title']?.[0] ||
      work.publisher ||
      'Unknown Journal';

    // Clean abstract (remove HTML/XML tags)
    const abstract = work.abstract
      ? work.abstract.replace(/<[^>]*>/g, '').trim()
      : '';

    // Extract DOI
    const doi = work.DOI || '';

    // Extract title
    const title = work.title?.[0] || 'Untitled';

    // Generate citation key
    const firstAuthorLastName = work.author?.[0]?.family || 'Unknown';
    const citationKey = `${firstAuthorLastName}_${year}`;

    return {
      title,
      authors,
      journal,
      year,
      doi,
      abstract,
      citationKey,
      itemType: work.type || 'journal-article',
      citationCount: work['is-referenced-by-count'],
      volume: work.volume,
      issue: work.issue,
      pages: work.page
    };
  } catch (error) {
    console.error('Error transforming work:', error);
    return null;
  }
}

/**
 * Generate BibTeX entry from publication data
 */
export function generateBibtex(pub: CrossrefPublication): string {
  const type = pub.itemType === 'journal-article' ? 'article' : 'misc';

  let bibtex = `@${type}{${pub.citationKey},\n`;
  bibtex += `  title = {${pub.title}},\n`;
  bibtex += `  author = {${pub.authors}},\n`;
  bibtex += `  journal = {${pub.journal}},\n`;
  bibtex += `  year = {${pub.year}},\n`;

  if (pub.volume) bibtex += `  volume = {${pub.volume}},\n`;
  if (pub.issue) bibtex += `  number = {${pub.issue}},\n`;
  if (pub.pages) bibtex += `  pages = {${pub.pages}},\n`;
  if (pub.doi) bibtex += `  doi = {${pub.doi}},\n`;

  bibtex += `}`;

  return bibtex;
}
