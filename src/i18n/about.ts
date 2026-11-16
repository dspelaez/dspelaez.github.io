import type { Lang } from './ui';

type Interest = { icon: string; title: string; description: string };
type Degree = { degree: string; institution: string; location: string; icon: string };

/** Short bio, in two boxes: what I do / my path. */
export const bio: Record<Lang, { what: string; path: string }> = {
  en: {
    what: 'I am a physical oceanographer specializing in air-sea interactions. I study ocean waves and how they influence the exchange of heat, momentum and gases in the upper layer of the ocean, combining field observations, satellite remote sensing and numerical models.',
    path: 'I trained in Colombia and Mexico, working on coastal extremes, coral reefs and sea turtle nesting beaches, then on air-sea turbulence. After a PhD in Ireland and France on different ways of measuring ocean waves, I now work as a postdoctoral researcher at LOPS in Brest, mapping ocean currents from space.',
  },
  es: {
    what: 'Soy oceanógrafo físico, especialista en interacciones océano-atmósfera. Estudio el oleaje y cómo influye en el intercambio de calor, momentum y gases en la capa superficial del océano, combinando observaciones de campo, teledetección satelital y modelos numéricos.',
    path: 'Me formé en Colombia y México, donde trabajé en eventos costeros extremos, arrecifes de coral y playas de anidación de tortugas, y luego en turbulencia océano-atmósfera. Tras un doctorado, en Irlanda y Francia, sobre diversas formas de medir el oleaje, hoy trabajo como investigador postdoctoral en el LOPS, en Brest, cartografiando corrientes oceánicas desde el espacio.',
  },
  fr: {
    what: 'Je suis océanographe physicien, spécialiste des interactions océan-atmosphère. J’étudie les vagues et leur influence sur les échanges de chaleur, de quantité de mouvement et de gaz dans la couche superficielle de l’océan, en combinant observations de terrain, télédétection satellitaire et modèles numériques.',
    path: 'Formé en Colombie et au Mexique, j’ai travaillé sur les événements côtiers extrêmes, les récifs coralliens et les plages de ponte des tortues, puis sur la turbulence océan-atmosphère. Après une thèse en Irlande et en France sur différentes manières de mesurer les vagues, je suis aujourd’hui chercheur postdoctoral au LOPS, à Brest, où je cartographie les courants océaniques depuis l’espace.',
  },
};

export const researchInterests: Record<Lang, Interest[]> = {
  en: [
    { icon: 'waves', title: 'Ocean wave dynamics', description: 'Generation, evolution, and dissipation processes' },
    { icon: 'wind', title: 'Air–sea momentum and energy transfer', description: 'Coupling between waves, wind, and currents' },
    { icon: 'turbulence', title: 'Upper ocean turbulence', description: 'Wave-driven boundary layers and Langmuir circulation' },
    { icon: 'bubbles', title: 'Gas exchange processes', description: 'Wave-induced and bubble-mediated air–sea fluxes' },
    { icon: 'sensors', title: 'Observational methods', description: 'Ocean acoustics, wave buoys, satellite remote sensing, and radar systems' },
    { icon: 'modeling', title: 'Numerical modeling', description: 'SWAN, WAVEWATCH III, Oceananigans.jl, GOTM, etc' },
  ],
  es: [
    { icon: 'waves', title: 'Dinámica del oleaje', description: 'Procesos de generación, evolución y disipación' },
    { icon: 'wind', title: 'Transferencia de momento y energía océano-atmósfera', description: 'Acoplamiento entre olas, viento y corrientes' },
    { icon: 'turbulence', title: 'Turbulencia en la capa superior del océano', description: 'Capas límite forzadas por el oleaje y circulación de Langmuir' },
    { icon: 'bubbles', title: 'Intercambio de gases', description: 'Flujos océano-atmósfera inducidos por el oleaje y mediados por burbujas' },
    { icon: 'sensors', title: 'Métodos de observación', description: 'Acústica oceánica, boyas de oleaje, teledetección satelital y sistemas de radar' },
    { icon: 'modeling', title: 'Modelado numérico', description: 'SWAN, WAVEWATCH III, Oceananigans.jl, GOTM, etc.' },
  ],
  fr: [
    { icon: 'waves', title: 'Dynamique des vagues', description: 'Processus de génération, d’évolution et de dissipation' },
    { icon: 'wind', title: 'Transferts de quantité de mouvement et d’énergie océan-atmosphère', description: 'Couplage entre vagues, vent et courants' },
    { icon: 'turbulence', title: 'Turbulence de la couche supérieure de l’océan', description: 'Couches limites forcées par les vagues et circulation de Langmuir' },
    { icon: 'bubbles', title: 'Échanges gazeux', description: 'Flux océan-atmosphère induits par les vagues et par les bulles' },
    { icon: 'sensors', title: 'Méthodes d’observation', description: 'Acoustique océanique, bouées houlographes, télédétection satellitaire et radars' },
    { icon: 'modeling', title: 'Modélisation numérique', description: 'SWAN, WAVEWATCH III, Oceananigans.jl, GOTM, etc.' },
  ],
};

export const education: Record<Lang, Degree[]> = {
  en: [
    { degree: 'PhD. Mathématiques aux interfaces', institution: 'Université Paris-Saclay', location: 'France', icon: 'phd' },
    { degree: 'PhD. Applied Mathematics', institution: 'University College Dublin', location: 'Ireland', icon: 'phd' },
    { degree: 'M.Sc. Physical Oceanography', institution: 'CICESE', location: 'Ensenada, México', icon: 'masters' },
    { degree: 'M.Eng. Water Resources', institution: 'Universidad Nacional de Colombia', location: 'Medellín', icon: 'masters' },
    { degree: 'B.Sc. Civil Engineering', institution: 'Universidad Nacional de Colombia', location: 'Medellín', icon: 'bachelor' },
  ],
  es: [
    { degree: 'Doctorado en Mathématiques aux interfaces', institution: 'Université Paris-Saclay', location: 'Francia', icon: 'phd' },
    { degree: 'Doctorado en Matemática Aplicada', institution: 'University College Dublin', location: 'Irlanda', icon: 'phd' },
    { degree: 'Maestría en Oceanografía Física', institution: 'CICESE', location: 'Ensenada, México', icon: 'masters' },
    { degree: 'Maestría en Ingeniería, Recursos Hidráulicos', institution: 'Universidad Nacional de Colombia', location: 'Medellín', icon: 'masters' },
    { degree: 'Ingeniería Civil', institution: 'Universidad Nacional de Colombia', location: 'Medellín', icon: 'bachelor' },
  ],
  fr: [
    { degree: 'Doctorat en mathématiques aux interfaces', institution: 'Université Paris-Saclay', location: 'France', icon: 'phd' },
    { degree: 'Doctorat en mathématiques appliquées', institution: 'University College Dublin', location: 'Irlande', icon: 'phd' },
    { degree: 'Master en océanographie physique', institution: 'CICESE', location: 'Ensenada, Mexique', icon: 'masters' },
    { degree: 'Master en ingénierie, ressources en eau', institution: 'Universidad Nacional de Colombia', location: 'Medellín, Colombie', icon: 'masters' },
    { degree: 'Diplôme d’ingénieur civil', institution: 'Universidad Nacional de Colombia', location: 'Medellín, Colombie', icon: 'bachelor' },
  ],
};
