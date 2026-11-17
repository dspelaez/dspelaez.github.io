/**
 * Blog posts and their translations.
 *
 * English posts live in src/content/blog/<name>.md and Spanish translations in
 * src/content/blog/es/<name>.md (same file name). A translation is linked to its
 * original through that shared name.
 */
import { getCollection, type CollectionEntry } from 'astro:content';
import type { Lang } from '../i18n/ui';

export type Post = CollectionEntry<'blog'>;

/** Languages that have translated posts (subfolder name). */
const translatedLangs = ['es', 'fr'] as const;

export function postLang(post: Post): Lang {
  const prefix = post.slug.split('/')[0];
  return (translatedLangs as readonly string[]).includes(prefix) ? (prefix as Lang) : 'en';
}

/** Post name without the language folder, shared by all its versions. */
export function postName(post: Post): string {
  return postLang(post) === 'en' ? post.slug : post.slug.split('/').slice(1).join('/');
}

export function postUrl(post: Post): string {
  const lang = postLang(post);
  return lang === 'en' ? `/blog/${postName(post)}/` : `/${lang}/blog/${postName(post)}/`;
}

/** All versions of a post, keyed by language (including the post itself). */
export async function postVersions(post: Post): Promise<Partial<Record<Lang, string>>> {
  const all = await getCollection('blog');
  const name = postName(post);
  return Object.fromEntries(
    all.filter((p) => postName(p) === name).map((p) => [postLang(p), postUrl(p)])
  );
}

export interface ListedPost {
  post: Post;
  url: string;
  /** Language the shown version is written in */
  shownLang: Lang;
  /** Shown in a language other than the page's (gets an EN/ES badge) */
  fallback: boolean;
  /** English post that also has a Spanish version */
  hasSpanish: boolean;
  /** The version picked for this page (one per article); the others are alternatives */
  primary: boolean;
  /** Topics of the article (set on the English original, shared by its translations) */
  topics: string[];
}

/**
 * Posts for a listing in the given language, newest first.
 *
 * Each article gets one primary entry: the version in the page's language when it
 * exists, otherwise the English one, otherwise any other (e.g. a post written only
 * in Spanish). With `allVersions`, the other language versions are included too
 * (primary: false), so a language filter can show every post written in a language.
 */
export async function listPosts(lang: Lang, { allVersions = false } = {}): Promise<ListedPost[]> {
  const all = await getCollection('blog');
  const articles = new Map<string, Partial<Record<Lang, Post>>>();
  for (const p of all) {
    const versions = articles.get(postName(p)) ?? {};
    versions[postLang(p)] = p;
    articles.set(postName(p), versions);
  }

  const entries: ListedPost[] = [];
  for (const versions of articles.values()) {
    // Page language first, then English, then Spanish, then French
    const preferred = (versions[lang] ?? versions.en ?? versions.es ?? versions.fr) as Post;
    const shown = allVersions ? (Object.values(versions) as Post[]) : [preferred];
    // Topics come from the original (English, or Spanish for posts written in Spanish)
    const topics = (versions.en ?? versions.es ?? preferred).data.topics;
    for (const post of shown) {
      const shownLang = postLang(post);
      entries.push({
        post,
        url: postUrl(post),
        shownLang,
        fallback: shownLang !== lang,
        hasSpanish: lang === 'en' && shownLang === 'en' && Boolean(versions.es),
        primary: post === preferred,
        topics,
      });
    }
  }
  return entries.sort((a, b) => b.post.data.pubDate.getTime() - a.post.data.pubDate.getTime());
}
