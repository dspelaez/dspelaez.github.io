---
title: "Una mejor forma de medir la dirección del oleaje: de la investigación al software"
description: "Cómo las wavelets nos ayudan a saber de dónde viene el oleaje, y por qué creamos una herramienta en Python para que cualquiera pueda usarlas."
pubDate: 2025-05-05
category: "research-update"
tags:
  - "Espectro direccional del oleaje"
  - "Wavelets"
  - "Boyas de oleaje"
  - "Python"
  - "Código abierto"
featured: true
---

## El problema de medir la dirección del oleaje

Imagina que estás en una playa durante un temporal. Las olas llegan desde varias direcciones: unas vienen de una borrasca lejana y otras las genera el viento local. Ahora imagina que intentas medir exactamente de dónde viene toda esa energía, frecuencia por frecuencia. Ese es el reto al que se enfrentan los oceanógrafos al analizar los datos de una boya de oleaje.

Durante décadas hemos usado métodos clásicos basados en el análisis de Fourier. Funcionan razonablemente bien, pero tienen limitaciones: suponen que el estado del mar no cambia durante la medida, a menudo exigen hacer hipótesis sobre la forma que debería tener el resultado y pueden fallar cuando se superponen varios sistemas de oleaje.

## Una solución basada en wavelets

En 2024 publicamos un artículo en el *Journal of Atmospheric and Oceanic Technology* ([Peláez-Zapata et al., 2024](https://doi.org/10.1175/JTECH-D-23-0058.1)) en el que proponíamos un enfoque basado en wavelets. Con datos de boyas GPS frente a la tempestuosa costa oeste de Irlanda, mostramos que las wavelets ofrecen ventajas reales:

**Localización en tiempo y frecuencia**: a diferencia de los métodos de Fourier, que analizan toda la serie temporal de una vez, las wavelets permiten fijarse a la vez en momentos y frecuencias concretas. Esto es clave cuando el estado del mar cambia rápidamente.

**Sin hipótesis previas**: el método se basa únicamente en los datos. No hace falta suponer que la distribución direccional tiene una forma determinada ni que las condiciones son estacionarias.

**Resultados robustos**: el método de wavelets produjo estimaciones direccionales más suaves y estables, incluso en temporales complejos con varios sistemas de oleaje superpuestos.

Validamos el método frente a los enfoques convencionales y a modelos numéricos, y comprobamos que funciona igual o mejor que los mejores métodos clásicos, sin sus hipótesis limitantes.

![Comparación de espectros direccionales del oleaje estimados con distintos métodos. El método de wavelets (panel a) produce distribuciones más suaves y realistas que los métodos clásicos (paneles b a e), y se parece mucho al modelo WaveWatch III (panel f).](/images/directional-spectrum-example.jpg)

*Figura: ejemplo de espectros direccionales del oleaje estimados con distintos métodos. El eje radial muestra la frecuencia (Hz) y la dirección angular indica hacia dónde se propaga el oleaje (N, E, S, W). Los colores indican la densidad de energía del oleaje para cada combinación de frecuencia y dirección.*

## Del artículo a la práctica: la herramienta EWDM

Los artículos científicos tienen un problema: demuestran que algo funciona, pero no siempre facilitan que otros lo usen. Por eso creamos **EWDM** (Extended Wavelet Directional Method), una herramienta de código abierto en Python que implementa estos algoritmos basados en wavelets.

EWDM se publicó en el *Journal of Open Source Software* ([Peláez-Zapata & Dias, 2025](https://doi.org/10.21105/joss.07942)) y pone esta investigación al alcance de toda la comunidad. La herramienta:

- Funciona con datos de fuentes muy diversas: boyas GPS, boyas de cabeceo y balanceo, sensores de nivel, ADCP e incluso vídeo estereoscópico.
- Se integra directamente con sistemas de boyas muy usados, como las boyas Spotter y la base de datos CDIP.
- Está construida sobre `xarray`, lo que la hace eficiente y compatible con los flujos de trabajo científicos modernos en Python.
- Incluye estimación de densidad por núcleos para evaluar de forma robusta la distribución direccional.

Está disponible en GitHub: [github.com/dspelaez/extended-wdm](https://github.com/dspelaez/extended-wdm)

## Por qué importa

Conocer bien la dirección del oleaje no es solo interesante desde el punto de vista académico; tiene aplicaciones muy concretas:

- **Seguridad marítima**: entender los mares cruzados y los patrones de oleaje complejos ayuda a evitar accidentes

- **Ingeniería de costas**: mejores datos direccionales permiten diseñar mejor puertos y defensas costeras

- **Energías renovables**: los parques eólicos marinos y los convertidores de energía del oleaje deben resistir olas que llegan desde varias direcciones

- **Investigación climática**: las estadísticas direccionales a largo plazo nos ayudan a entender cómo está cambiando el oleaje

## Lo que aprendimos en la costa irlandesa

La costa oeste de Irlanda fue el lugar ideal para poner a prueba el método. Esta región recibe toda la fuerza de los temporales del Atlántico Norte, con un oleaje complejo que cambia rápidamente. Un resultado especialmente interesante fue que la dispersión direccional (cuánto se reparten las olas alrededor de la dirección principal) apenas dependía de la edad del oleaje. Esto sugiere que las interacciones no lineales entre olas influyen más que el forzamiento directo del viento en cómo se reparte el oleaje en dirección. También encontramos que la dispersión direccional era mínima en el pico espectral (las olas dominantes) y aumentaba de forma asimétrica hacia frecuencias más altas y más bajas.

## Ciencia abierta en acción

Este proyecto encarna lo que debería ser la ciencia abierta: publicar la investigación, liberar el código y facilitar que otros construyan sobre tu trabajo. EWDM tiene una licencia de código abierto permisiva, está archivado en Zenodo y cuenta con documentación y ejemplos.

Tanto si procesas los datos de una sola boya como si analizas series largas de clima de oleaje, EWDM ofrece una forma flexible y sin hipótesis previas de estimar el espectro direccional.

## Para empezar

Si trabajas con datos de oleaje y quieres probar el método de wavelets:

```bash
pip install ewdm
```

El [repositorio en GitHub](https://github.com/dspelaez/extended-wdm) incluye documentación, ejemplos y notebooks de Jupyter para empezar rápidamente.

## Reflexión final

A veces el mejor avance científico surge de combinar herramientas matemáticas que ya existían (las wavelets llevan décadas entre nosotros) con el conocimiento del dominio (entender el oleaje) y el compromiso de hacer accesibles los resultados. Eso es lo que representa este proyecto.

Y si alguna vez estás en la costa irlandesa viendo llegar un mar de fondo enorme durante un temporal de invierno, podrás apreciar la sofisticación matemática que hace falta para medir de dónde viene realmente toda esa energía, y saber que ahora existe una herramienta gratuita para hacerlo.

---

**Referencias:**

- Peláez-Zapata, D.S., et al. (2024). Ocean Wave Directional Distribution from GPS Buoy Observations off the West Coast of Ireland: Assessment of a Wavelet-Based Method. *J. Atmos. Ocean. Technol.* [doi: 10.1175/JTECH-D-23-0058.1](https://doi.org/10.1175/JTECH-D-23-0058.1)

- Peláez-Zapata, D.S., & Dias, F. (2025). EWDM: A wavelet-based method for estimating directional spectra of ocean waves. *J. Open Source Softw.*, 10(109), 7942. [doi: 10.21105/joss.07942](https://doi.org/10.21105/joss.07942)
