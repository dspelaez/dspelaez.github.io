export const languages = {
  en: 'English',
  es: 'Español',
  fr: 'Français',
} as const;

export type Lang = keyof typeof languages;
export const defaultLang: Lang = 'en';
export const langs = Object.keys(languages) as Lang[];

/** Locale used for dates and Open Graph. */
export const locales: Record<Lang, string> = {
  en: 'en-US',
  es: 'es-ES',
  fr: 'fr-FR',
};

/** Language of the current page, from its URL (/es/..., /fr/..., otherwise English). */
export function getLang(url: URL): Lang {
  const first = url.pathname.split('/')[1];
  return (langs as string[]).includes(first) ? (first as Lang) : defaultLang;
}

/** Path without its language prefix, e.g. /es/blog/ -> /blog/ */
export function stripLang(pathname: string): string {
  const parts = pathname.split('/');
  if ((langs as string[]).includes(parts[1])) parts.splice(1, 1);
  return parts.join('/') || '/';
}

/** Path for a given language, e.g. ('/blog/', 'fr') -> /fr/blog/ */
export function localizePath(path: string, lang: Lang): string {
  const clean = stripLang(path);
  return lang === defaultLang ? clean : `/${lang}${clean === '/' ? '/' : clean}`;
}

/** Pages that exist in every language (paths without prefix). */
export const translatedPages = ['/', '/blog/', '/projects/', '/publications/'];

export function isTranslatedPage(pathname: string): boolean {
  const base = stripLang(pathname).replace(/\/?$/, '/');
  return translatedPages.includes(base);
}

export function formatDate(date: Date, lang: Lang, opts: Intl.DateTimeFormatOptions) {
  return date.toLocaleDateString(locales[lang], opts);
}

/** French typography puts a space before the colon. */
export const colon = (lang: Lang) => (lang === 'fr' ? ' :' : ':');

const ui = {
  en: {
    'meta.home.title': 'Daniel Peláez-Zapata | Physical Oceanographer | PhD Applied Mathematics',
    'meta.home.description': 'Physical oceanographer specializing in ocean wave dynamics, air-sea interactions, and numerical modeling. Research on energy transfer, momentum exchange, and mass fluxes between ocean and atmosphere.',
    'meta.blog.title': 'Blog & Tutorials | Daniel Peláez-Zapata',
    'meta.blog.description': 'Technical guides, research updates, and tutorials on Python, SWAN, WAVEWATCH III, and data analysis',
    'meta.projects.title': 'Research Projects | Daniel Peláez-Zapata',
    'meta.projects.description': 'Current and past research projects in physical oceanography, wave modeling, and coastal dynamics',
    'meta.publications.title': 'Publications | Daniel Peláez-Zapata',
    'meta.publications.description': 'Complete list of peer-reviewed publications in physical oceanography, ocean waves, and coastal dynamics',

    'nav.about': 'About',
    'nav.publications': 'Publications',
    'nav.projects': 'Projects',
    'nav.blog': 'Blog',
    'nav.contact': 'Contact',
    'nav.theme': 'Switch between light and dark mode',
    'nav.language': 'Language',
    'nav.menu': 'Toggle menu',

    'hero.subtitle': 'Researcher in physical oceanography and air-sea interactions',
    'hero.position': 'Postdoctoral researcher · LOPS (Ifremer / CNRS)',

    'about.title': 'About Me',
    'about.what': 'What I do',
    'about.path': 'My path',
    'about.interests': 'Research Interests',
    'about.education': 'Education',

    'carousel.prev': 'Previous',
    'carousel.next': 'Next',
    'carousel.goto': 'Go to slide',

    'pubs.title': 'Recent Publications',
    'pubs.subtitle': 'Selected peer-reviewed publications in physical oceanography',
    'pubs.pageTitle': 'Publications',
    'pubs.pageSubtitle': 'Peer-reviewed publications in oceanography journals',
    'pubs.viewAll': 'View All Publications',
    'pubs.citation': 'citation',
    'pubs.citations': 'citations',
    'pubs.link': 'Link',
    'pubs.copy': 'Copy BibTeX',
    'pubs.copied': '✓ Copied!',
    'pubs.copyFailed': 'Failed to copy BibTeX to clipboard',
    'pubs.showAbstract': 'Show abstract',
    'pubs.conferenceAbstract': 'Conference abstract',
    'pubs.phdThesis': 'PhD thesis',

    'projects.title': 'Research Projects',
    'projects.subtitle': 'Current and future research initiatives',
    'projects.pageSubtitle': 'Ongoing, planned, and completed research initiatives',
    'projects.ongoing': 'Ongoing Projects',
    'projects.planned': 'Planned Projects',
    'projects.completed': 'Completed Projects',
    'projects.status.ongoing': 'ongoing',
    'projects.status.planned': 'planned',
    'projects.status.completed': 'completed',
    'projects.funding': 'Funding',
    'projects.pi': 'PI',
    'projects.pis': 'PIs',
    'projects.collaborators': 'Collaborators',
    'projects.started': 'Started',
    'projects.start': 'Start',
    'projects.repo': 'View Repository',
    'projects.viewAll': 'View All Projects',

    'blog.title': 'Blog & Tutorials',
    'blog.subtitle': 'Technical guides, research updates, and insights',
    'blog.pageSubtitle': 'Technical guides, research updates, and oceanography tutorials',
    'blog.readMore': 'Read more',
    'blog.back': 'Back to Blog',
    'blog.backAll': 'Back to all posts',
    'blog.readIn': 'Read this post in English',
    'blog.alsoSpanish': 'También en español',
    'blog.only.en': 'Available in English only',
    'blog.only.es': 'Available in Spanish only',
    'blog.viewAll': 'View All Posts',
    'blog.englishOnly': '',
    'blog.cat.research-update': 'research update',
    'blog.cat.tutorial': 'tutorial',
    'blog.filter.all': 'All',
    'blog.filter.research-update': 'Research updates',
    'blog.filter.tutorial': 'Tutorials',
    'blog.filter.lang.all': 'All',
    'blog.filter.lang.en': 'English',
    'blog.filter.lang.es': 'Spanish',
    'blog.filter.langLabel': 'Language',
    'blog.filter.topicLabel': 'Topics',
    'blog.topic.wave-measurements': 'Wave measurements',
    'blog.topic.air-sea-climate': 'Air-sea interaction & climate',
    'blog.topic.coastal-hazards': 'Coastal hazards',
    'blog.topic.marine-energy': 'Marine renewable energy',
    'blog.topic.modelling': 'Modelling & computing',
    'blog.filter.label': 'Filter posts',
    'blog.filter.empty': 'No posts match this filter.',

    'contact.title': 'Get in Touch',
    'contact.subtitle': 'Connect with me on various platforms',
    'contact.profile': 'Profile',
    'contact.collabTitle': 'Collaboration & Opportunities',
    'contact.collabText': "I'm always interested in collaborating on research projects, consulting opportunities, or discussing ocean science. Feel free to reach out!",
    'contact.send': 'Send Email',

    'footer.rights': 'All rights reserved.',
    'footer.builtWith': 'Built with',
    'footer.and': 'and',
  },
  es: {
    'meta.home.title': 'Daniel Peláez-Zapata | Oceanógrafo físico | Doctor en Matemática Aplicada',
    'meta.home.description': 'Oceanógrafo físico especializado en dinámica del oleaje, interacción océano-atmósfera y modelado numérico. Investigación sobre la transferencia de energía, momento y masa entre el océano y la atmósfera.',
    'meta.blog.title': 'Blog y tutoriales | Daniel Peláez-Zapata',
    'meta.blog.description': 'Guías técnicas, novedades de investigación y tutoriales sobre Python, SWAN, WAVEWATCH III y análisis de datos',
    'meta.projects.title': 'Proyectos de investigación | Daniel Peláez-Zapata',
    'meta.projects.description': 'Proyectos de investigación actuales y pasados en oceanografía física, modelado del oleaje y dinámica costera',
    'meta.publications.title': 'Publicaciones | Daniel Peláez-Zapata',
    'meta.publications.description': 'Lista completa de publicaciones arbitradas en oceanografía física, oleaje y dinámica costera',

    'nav.about': 'Sobre mí',
    'nav.publications': 'Publicaciones',
    'nav.projects': 'Proyectos',
    'nav.blog': 'Blog',
    'nav.contact': 'Contacto',
    'nav.theme': 'Cambiar entre modo claro y oscuro',
    'nav.language': 'Idioma',
    'nav.menu': 'Abrir o cerrar el menú',

    'hero.subtitle': 'Investigador en oceanografía física e interacción océano-atmósfera',
    'hero.position': 'Investigador posdoctoral · LOPS (Ifremer / CNRS)',

    'about.title': 'Sobre mí',
    'about.what': '¿A qué me dedico?',
    'about.path': 'Mi recorrido',
    'about.interests': 'Intereses de investigación',
    'about.education': 'Formación',

    'carousel.prev': 'Anterior',
    'carousel.next': 'Siguiente',
    'carousel.goto': 'Ir a la diapositiva',

    'pubs.title': 'Publicaciones recientes',
    'pubs.subtitle': 'Selección de publicaciones arbitradas en oceanografía física',
    'pubs.pageTitle': 'Publicaciones',
    'pubs.pageSubtitle': 'Publicaciones arbitradas en revistas de oceanografía',
    'pubs.viewAll': 'Ver todas las publicaciones',
    'pubs.citation': 'cita',
    'pubs.citations': 'citas',
    'pubs.link': 'Enlace',
    'pubs.copy': 'Copiar BibTeX',
    'pubs.copied': '✓ ¡Copiado!',
    'pubs.copyFailed': 'No se pudo copiar el BibTeX',
    'pubs.showAbstract': 'Ver resumen',
    'pubs.conferenceAbstract': 'Resumen de congreso',
    'pubs.phdThesis': 'Tesis doctoral',

    'projects.title': 'Proyectos de investigación',
    'projects.subtitle': 'Iniciativas de investigación actuales y futuras',
    'projects.pageSubtitle': 'Iniciativas de investigación en curso, planificadas y finalizadas',
    'projects.ongoing': 'Proyectos en curso',
    'projects.planned': 'Proyectos planificados',
    'projects.completed': 'Proyectos finalizados',
    'projects.status.ongoing': 'en curso',
    'projects.status.planned': 'planificado',
    'projects.status.completed': 'finalizado',
    'projects.funding': 'Financiación',
    'projects.pi': 'Investigador principal',
    'projects.pis': 'Investigadores principales',
    'projects.collaborators': 'Colaboradores',
    'projects.started': 'Inicio',
    'projects.start': 'Inicio',
    'projects.repo': 'Ver repositorio',
    'projects.viewAll': 'Ver todos los proyectos',

    'blog.title': 'Blog y tutoriales',
    'blog.subtitle': 'Guías técnicas, novedades de investigación y reflexiones',
    'blog.pageSubtitle': 'Guías técnicas, novedades de investigación y tutoriales de oceanografía',
    'blog.readMore': 'Leer más',
    'blog.back': 'Volver al blog',
    'blog.backAll': 'Volver a todos los artículos',
    'blog.readIn': 'Lee este artículo en español',
    'blog.alsoSpanish': 'También en español',
    'blog.only.en': 'Disponible solo en inglés',
    'blog.only.es': 'Disponible solo en español',
    'blog.viewAll': 'Ver todos los artículos',
    'blog.englishOnly': 'Algunos artículos están disponibles solo en inglés (marcados con EN).',
    'blog.cat.research-update': 'novedades',
    'blog.cat.tutorial': 'tutorial',
    'blog.filter.all': 'Todos',
    'blog.filter.research-update': 'Novedades',
    'blog.filter.tutorial': 'Tutoriales',
    'blog.filter.lang.all': 'Todos',
    'blog.filter.lang.en': 'Inglés',
    'blog.filter.lang.es': 'Español',
    'blog.filter.langLabel': 'Idioma',
    'blog.filter.topicLabel': 'Temas',
    'blog.topic.wave-measurements': 'Medición del oleaje',
    'blog.topic.air-sea-climate': 'Océano-atmósfera y clima',
    'blog.topic.coastal-hazards': 'Riesgos costeros',
    'blog.topic.marine-energy': 'Energías marinas',
    'blog.topic.modelling': 'Modelado y computación',
    'blog.filter.label': 'Filtrar artículos',
    'blog.filter.empty': 'Ningún artículo coincide con este filtro.',

    'contact.title': 'Contacto',
    'contact.subtitle': 'Encuéntrame en estas plataformas',
    'contact.profile': 'Perfil',
    'contact.collabTitle': 'Colaboraciones y oportunidades',
    'contact.collabText': 'Siempre me interesa colaborar en proyectos de investigación, oportunidades de consultoría o conversar sobre ciencias del mar. ¡No dudes en escribirme!',
    'contact.send': 'Enviar correo',

    'footer.rights': 'Todos los derechos reservados.',
    'footer.builtWith': 'Hecho con',
    'footer.and': 'y',
  },
  fr: {
    'meta.home.title': 'Daniel Peláez-Zapata | Océanographe physicien | Docteur en mathématiques appliquées',
    'meta.home.description': "Océanographe physicien spécialisé dans la dynamique des vagues, les interactions océan-atmosphère et la modélisation numérique. Recherches sur les transferts d'énergie, de quantité de mouvement et de masse entre l'océan et l'atmosphère.",
    'meta.blog.title': 'Blog et tutoriels | Daniel Peláez-Zapata',
    'meta.blog.description': 'Guides techniques, actualités de la recherche et tutoriels sur Python, SWAN, WAVEWATCH III et l’analyse de données',
    'meta.projects.title': 'Projets de recherche | Daniel Peláez-Zapata',
    'meta.projects.description': 'Projets de recherche actuels et passés en océanographie physique, modélisation des vagues et dynamique côtière',
    'meta.publications.title': 'Publications | Daniel Peláez-Zapata',
    'meta.publications.description': 'Liste complète des publications à comité de lecture en océanographie physique, vagues et dynamique côtière',

    'nav.about': 'À propos',
    'nav.publications': 'Publications',
    'nav.projects': 'Projets',
    'nav.blog': 'Blog',
    'nav.contact': 'Contact',
    'nav.theme': 'Basculer entre mode clair et mode sombre',
    'nav.language': 'Langue',
    'nav.menu': 'Ouvrir ou fermer le menu',

    'hero.subtitle': 'Chercheur en océanographie physique et interactions océan-atmosphère',
    'hero.position': 'Chercheur postdoctoral · LOPS (Ifremer / CNRS)',

    'about.title': 'À propos',
    'about.what': 'Ce que je fais',
    'about.path': 'Mon parcours',
    'about.interests': 'Thèmes de recherche',
    'about.education': 'Formation',

    'carousel.prev': 'Précédent',
    'carousel.next': 'Suivant',
    'carousel.goto': 'Aller à la diapositive',

    'pubs.title': 'Publications récentes',
    'pubs.subtitle': 'Sélection de publications à comité de lecture en océanographie physique',
    'pubs.pageTitle': 'Publications',
    'pubs.pageSubtitle': "Publications à comité de lecture dans des revues d'océanographie",
    'pubs.viewAll': 'Voir toutes les publications',
    'pubs.citation': 'citation',
    'pubs.citations': 'citations',
    'pubs.link': 'Lien',
    'pubs.copy': 'Copier le BibTeX',
    'pubs.copied': '✓ Copié !',
    'pubs.copyFailed': 'Impossible de copier le BibTeX',
    'pubs.showAbstract': 'Voir le résumé',
    'pubs.conferenceAbstract': 'Résumé de conférence',
    'pubs.phdThesis': 'Thèse de doctorat',

    'projects.title': 'Projets de recherche',
    'projects.subtitle': 'Initiatives de recherche actuelles et futures',
    'projects.pageSubtitle': 'Initiatives de recherche en cours, prévues et terminées',
    'projects.ongoing': 'Projets en cours',
    'projects.planned': 'Projets prévus',
    'projects.completed': 'Projets terminés',
    'projects.status.ongoing': 'en cours',
    'projects.status.planned': 'prévu',
    'projects.status.completed': 'terminé',
    'projects.funding': 'Financement',
    'projects.pi': 'Responsable scientifique',
    'projects.pis': 'Responsables scientifiques',
    'projects.collaborators': 'Collaborateurs',
    'projects.started': 'Début',
    'projects.start': 'Début',
    'projects.repo': 'Voir le dépôt',
    'projects.viewAll': 'Voir tous les projets',

    'blog.title': 'Blog et tutoriels',
    'blog.subtitle': 'Guides techniques, actualités de la recherche et réflexions',
    'blog.pageSubtitle': "Guides techniques, actualités de la recherche et tutoriels d'océanographie",
    'blog.readMore': 'Lire la suite',
    'blog.back': 'Retour au blog',
    'blog.backAll': 'Retour à tous les articles',
    'blog.readIn': 'Lire cet article en français',
    'blog.alsoSpanish': 'También en español',
    'blog.only.en': 'Disponible uniquement en anglais',
    'blog.only.es': 'Disponible uniquement en espagnol',
    'blog.viewAll': 'Voir tous les articles',
    'blog.englishOnly': "Les articles ne sont pas encore traduits en français (langue indiquée par EN ou ES).",
    'blog.cat.research-update': 'actualité',
    'blog.cat.tutorial': 'tutoriel',
    'blog.filter.all': 'Tous',
    'blog.filter.research-update': 'Actualités',
    'blog.filter.tutorial': 'Tutoriels',
    'blog.filter.lang.all': 'Toutes',
    'blog.filter.lang.en': 'Anglais',
    'blog.filter.lang.es': 'Espagnol',
    'blog.filter.langLabel': 'Langue',
    'blog.filter.topicLabel': 'Thèmes',
    'blog.topic.wave-measurements': 'Mesure des vagues',
    'blog.topic.air-sea-climate': 'Océan-atmosphère et climat',
    'blog.topic.coastal-hazards': 'Risques côtiers',
    'blog.topic.marine-energy': 'Énergies marines',
    'blog.topic.modelling': 'Modélisation et calcul',
    'blog.filter.label': 'Filtrer les articles',
    'blog.filter.empty': 'Aucun article ne correspond à ce filtre.',

    'contact.title': 'Contact',
    'contact.subtitle': 'Retrouvez-moi sur ces plateformes',
    'contact.profile': 'Profil',
    'contact.collabTitle': 'Collaborations et opportunités',
    'contact.collabText': "Je suis toujours intéressé par des collaborations de recherche, des missions de conseil ou des échanges sur les sciences de la mer. N'hésitez pas à me contacter !",
    'contact.send': 'Envoyer un e-mail',

    'footer.rights': 'Tous droits réservés.',
    'footer.builtWith': 'Réalisé avec',
    'footer.and': 'et',
  },
} as const;

export type UIKey = keyof (typeof ui)['en'];

export function useTranslations(lang: Lang) {
  // Accepts built keys like `projects.status.${status}`; falls back to English
  return (key: UIKey | string): string => {
    const k = key as UIKey;
    return (ui[lang] as Record<UIKey, string>)[k] ?? ui[defaultLang][k] ?? key;
  };
}

/** Venue labels that come in English from the publication metadata. */
export function venueLabel(journal: string, lang: Lang): string {
  const t = useTranslations(lang);
  if (journal === 'Conference abstract') return t('pubs.conferenceAbstract');
  if (journal.startsWith('PhD thesis')) return journal.replace('PhD thesis', t('pubs.phdThesis'));
  return journal;
}

/** Title and description of a project in the given language (English if missing). */
export function projectText(
  project: { data: { title: string; description: string; i18n?: Partial<Record<Lang, { title?: string; description?: string }>> } },
  lang: Lang
) {
  const tr = project.data.i18n?.[lang];
  return {
    title: tr?.title ?? project.data.title,
    description: tr?.description ?? project.data.description,
  };
}
