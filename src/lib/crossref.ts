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
  url?: string;
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
  relation?: Record<string, unknown[]>;
  institution?: { name?: string }[];
  approved?: { 'date-parts'?: number[][] };
}

interface OrcidWorkSummary {
  'put-code': number;
  type?: string;
  title?: { title?: { value?: string } };
  'journal-title'?: { value?: string } | null;
  'publication-date'?: { year?: { value?: string } | null } | null;
  'external-ids'?: { 'external-id'?: { 'external-id-type': string; 'external-id-value': string }[] };
  url?: { value?: string } | null;
}

const ORCID_API = 'https://pub.orcid.org/v3.0';

const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

async function getJson(url: string, headers: Record<string, string> = {}, retries = 3) {
  const response = await fetch(url, { headers });
  // Rate limited: wait and try again
  if (response.status === 429 && retries > 0) {
    await sleep(2000);
    return getJson(url, headers, retries - 1);
  }
  if (!response.ok) {
    throw new Error(`${url}: ${response.status} ${response.statusText}`);
  }
  return response.json();
}

/**
 * Fetch publications listed in an ORCID record.
 *
 * ORCID gives the list of works (the user controls it), and CrossRef adds the
 * details (authors, citations, volume) for works that have a DOI. Asking
 * CrossRef directly for an ORCID misses older papers whose metadata never
 * recorded the author's ORCID.
 */
export async function fetchPublications(
  orcid: string,
  email: string = 'noreply@example.com'
): Promise<CrossrefPublication[]> {
  console.log(`Fetching works from ORCID: ${orcid}`);
  const record = await getJson(`${ORCID_API}/${orcid}/works`, { Accept: 'application/json' });
  const summaries: OrcidWorkSummary[] = (record.group || []).map(
    (g: { 'work-summary': OrcidWorkSummary[] }) => g['work-summary'][0]
  );
  console.log(`Found ${summaries.length} works in ORCID`);

  // One request at a time: CrossRef rate-limits bursts of parallel requests
  const results: (CrossrefPublication | null)[] = [];
  for (const summary of summaries) {
    results.push(await (async () => {
      const doi = summary['external-ids']?.['external-id']?.find(
        (e) => e['external-id-type'] === 'doi'
      )?.['external-id-value'];
      const orcidYear = Number(summary['publication-date']?.year?.value) || undefined;

      if (doi) {
        try {
          const data = await getJson(
            `https://api.crossref.org/works/${encodeURIComponent(doi)}?mailto=${email}`
          );
          const work: CrossrefWork = data.message;
          // Skip preprints that already have a published version
          if (work.relation?.['is-preprint-of']?.length) return null;
          const pub = transformWork(work);
          if (pub) {
            // The ORCID record is the reference for the year
            if (orcidYear) pub.year = orcidYear;
            pub.url = `https://doi.org/${pub.doi}`;
          }
          return pub;
        } catch (error) {
          console.warn(`CrossRef lookup failed for ${doi}, using ORCID data`, error);
        }
      }

      return orcidWork(orcid, summary, doi, orcidYear);
    })());
  }

  return (results.filter(Boolean) as CrossrefPublication[]).sort((a, b) => b.year - a.year);
}

/**
 * Build a publication from ORCID alone (works without DOI, or CrossRef down).
 */
async function orcidWork(
  orcid: string,
  summary: OrcidWorkSummary,
  doi: string | undefined,
  year: number | undefined
): Promise<CrossrefPublication | null> {
  const title = summary.title?.title?.value;
  if (!title) return null;

  let authors = '';
  try {
    const detail = await getJson(`${ORCID_API}/${orcid}/work/${summary['put-code']}`, {
      Accept: 'application/json',
    });
    authors = (detail.contributors?.contributor || [])
      .map((c: { 'credit-name'?: { value?: string } }) => c['credit-name']?.value)
      .filter(Boolean)
      .join(', ');
  } catch (error) {
    console.warn(`ORCID detail failed for "${title}"`, error);
  }

  const finalYear = year || new Date().getFullYear();
  const firstFamily = authors.split(',')[0]?.trim().split(' ').pop() || 'Unknown';

  return {
    title,
    authors: authors || 'Unknown',
    journal: summary['journal-title']?.value || '',
    year: finalYear,
    doi: doi || '',
    abstract: '',
    citationKey: `${firstFamily}_${finalYear}`,
    itemType: summary.type || 'other',
    url: doi ? `https://doi.org/${doi}` : summary.url?.value,
  };
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
      work.approved?.['date-parts']?.[0]?.[0] ||
      new Date().getFullYear();

    // Theses: show the degree-granting institution instead of the registrar
    const institution = work.institution?.[0]?.name;
    const thesisVenue = work.type === 'dissertation'
      ? `PhD thesis${institution ? `, ${institution.charAt(0).toUpperCase()}${institution.slice(1)}` : ''}`
      : '';

    // Extract journal name (posted-content left after the preprint filter
    // are conference abstracts, e.g. EGU/EMS, published by Copernicus)
    const journal =
      (work.type === 'posted-content' ? 'Conference abstract' : '') ||
      thesisVenue ||
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
  const type =
    pub.itemType === 'journal-article' ? 'article'
    : pub.itemType === 'conference-paper' ? 'inproceedings'
    : pub.itemType === 'dissertation' ? 'phdthesis'
    : 'misc';
  const venueField = type === 'inproceedings' ? 'booktitle' : type === 'phdthesis' ? 'school' : 'journal';

  let bibtex = `@${type}{${pub.citationKey},\n`;
  bibtex += `  title = {${pub.title}},\n`;
  bibtex += `  author = {${pub.authors}},\n`;
  if (pub.journal) bibtex += `  ${venueField} = {${pub.journal.replace(/^PhD thesis, /, '')}},\n`;
  bibtex += `  year = {${pub.year}},\n`;

  if (pub.volume) bibtex += `  volume = {${pub.volume}},\n`;
  if (pub.issue) bibtex += `  number = {${pub.issue}},\n`;
  if (pub.pages) bibtex += `  pages = {${pub.pages}},\n`;
  if (pub.doi) bibtex += `  doi = {${pub.doi}},\n`;
  else if (pub.url) bibtex += `  url = {${pub.url}},\n`;

  bibtex += `}`;

  return bibtex;
}
