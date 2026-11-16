/**
 * Structured data (schema.org JSON-LD) so search engines and AI assistants can
 * identify the site's author as a single person across name variants and profiles.
 */
import type { Lang } from '../i18n/ui';

const SITE = 'https://dspelaez.github.io';
export const PERSON_ID = `${SITE}/#person`;

const description: Record<Lang, string> = {
  en: 'Physical oceanographer specializing in air-sea interactions and ocean waves. Postdoctoral researcher at LOPS (Ifremer / CNRS), Brest.',
  es: 'Oceanógrafo físico especialista en interacciones océano-atmósfera y oleaje. Investigador posdoctoral en el LOPS (Ifremer / CNRS), Brest.',
  fr: 'Océanographe physicien spécialiste des interactions océan-atmosphère et des vagues. Chercheur postdoctoral au LOPS (Ifremer / CNRS), Brest.',
};

const org = (name: string, url?: string) => ({ '@type': 'Organization', name, ...(url && { url }) });

export function personSchema(lang: Lang) {
  return {
    '@context': 'https://schema.org',
    '@type': 'Person',
    '@id': PERSON_ID,
    name: 'Daniel Peláez-Zapata',
    // Spellings found in publication metadata
    alternateName: ['Daniel Santiago Peláez Zapata', 'Daniel Pelaez-Zapata', 'Daniel Pelaez Zapata', 'D. S. Peláez-Zapata'],
    url: `${SITE}/`,
    image: `${SITE}/profile.jpg`,
    email: 'mailto:dspelaez@ifremer.fr',
    jobTitle: lang === 'es' ? 'Investigador posdoctoral' : lang === 'fr' ? 'Chercheur postdoctoral' : 'Postdoctoral researcher',
    description: description[lang],
    worksFor: {
      '@type': 'Organization',
      name: "Laboratoire d'Océanographie Physique et Spatiale (LOPS)",
      parentOrganization: [org('Ifremer', 'https://www.ifremer.fr'), org('CNRS', 'https://www.cnrs.fr')],
      address: { '@type': 'PostalAddress', addressLocality: 'Brest', addressCountry: 'FR' },
    },
    alumniOf: [
      org('Université Paris-Saclay'),
      org('University College Dublin'),
      org('CICESE'),
      org('Universidad Nacional de Colombia'),
    ],
    knowsAbout: [
      'Physical oceanography',
      'Ocean waves',
      'Air-sea interactions',
      'Upper ocean turbulence',
      'Satellite remote sensing',
      'Directional wave spectra',
      'Meteotsunamis',
      'Numerical modelling',
    ],
    knowsLanguage: ['es', 'en', 'fr'],
    nationality: { '@type': 'Country', name: 'Colombia' },
    sameAs: [
      'https://orcid.org/0000-0001-5862-6194',
      'https://scholar.google.com/citations?user=Usm9u2AAAAAJ',
      'https://github.com/dspelaez',
      'https://www.researchgate.net/profile/Daniel-Pelaez-Zapata-3',
    ],
  };
}

export function websiteSchema(lang: Lang) {
  return {
    '@context': 'https://schema.org',
    '@type': 'WebSite',
    '@id': `${SITE}/#website`,
    url: `${SITE}/`,
    name: 'Daniel Peláez-Zapata',
    inLanguage: ['en', 'es', 'fr'],
    author: { '@id': PERSON_ID },
    description: description[lang],
  };
}

export function blogPostingSchema(post: {
  title: string;
  description: string;
  url: string;
  pubDate: Date;
  image?: string;
  tags?: string[];
  paperUrl?: string;
}) {
  return {
    '@context': 'https://schema.org',
    '@type': 'BlogPosting',
    headline: post.title,
    description: post.description,
    url: post.url,
    mainEntityOfPage: post.url,
    datePublished: post.pubDate.toISOString(),
    inLanguage: 'en',
    author: { '@type': 'Person', '@id': PERSON_ID, name: 'Daniel Peláez-Zapata', url: `${SITE}/` },
    ...(post.image && { image: new URL(post.image, SITE).href }),
    ...(post.tags?.length && { keywords: post.tags.join(', ') }),
    // Posts written about one of the author's papers
    ...(post.paperUrl && { isBasedOn: { '@type': 'ScholarlyArticle', url: post.paperUrl } }),
  };
}
