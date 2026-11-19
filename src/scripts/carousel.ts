/**
 * Homepage carousels (blog and projects).
 *
 * Expects #<name>-carousel (the scrolling container), .<name>-carousel-item cards,
 * .<name>-carousel-indicator dots (one per card in the markup) and
 * .<name>-carousel-prev / -next buttons.
 *
 * There is one position per possible scroll step, not per card: with 6 cards and 3
 * visible there are 4 positions (cards 1-3, 2-4, 3-5, 4-6). Only that many dots are
 * shown, and the carousel goes back to the start after the last position. The
 * count is recomputed when the window is resized (3, 2 or 1 cards visible).
 */
const AUTO_SCROLL_DELAY = 5000; // ms

export function initCarousel(name: string) {
  const carousel = document.getElementById(`${name}-carousel`);
  const items = document.querySelectorAll<HTMLElement>(`.${name}-carousel-item`);
  const indicators = document.querySelectorAll<HTMLButtonElement>(`.${name}-carousel-indicator`);
  const prevBtn = document.querySelector<HTMLButtonElement>(`.${name}-carousel-prev`);
  const nextBtn = document.querySelector<HTMLButtonElement>(`.${name}-carousel-next`);

  if (!carousel || items.length === 0) return;

  let currentIndex = 0;
  let lastIndex = 0;
  let autoScroll: ReturnType<typeof setInterval> | undefined;

  // Distance between the start of two consecutive cards
  const step = () => {
    const gap = parseFloat(getComputedStyle(items[0].parentElement!).columnGap) || 24;
    return items[0].offsetWidth + gap;
  };

  // Last position: how many steps the container can actually scroll
  const computeLastIndex = () =>
    Math.max(0, Math.round((carousel.scrollWidth - carousel.clientWidth) / step()));

  function updateIndicators() {
    indicators.forEach((dot, i) => {
      dot.hidden = i > lastIndex;
      const active = i === currentIndex;
      dot.classList.toggle('bg-blue-600', active);
      dot.classList.toggle('dark:bg-blue-500', active);
      dot.classList.toggle('w-8', active);
      dot.classList.toggle('bg-slate-300', !active);
      dot.classList.toggle('dark:bg-slate-600', !active);
      dot.setAttribute('aria-current', active ? 'true' : 'false');
    });
  }

  function goTo(index: number, behavior: ScrollBehavior = 'smooth') {
    currentIndex = Math.min(Math.max(index, 0), lastIndex);
    carousel!.scrollTo({ left: step() * currentIndex, behavior });
    updateIndicators();
  }

  const next = () => goTo(currentIndex >= lastIndex ? 0 : currentIndex + 1);
  const prev = () => goTo(currentIndex <= 0 ? lastIndex : currentIndex - 1);

  function stopAutoScroll() {
    clearInterval(autoScroll);
    autoScroll = undefined;
  }
  function startAutoScroll() {
    stopAutoScroll();
    autoScroll = setInterval(next, AUTO_SCROLL_DELAY);
  }
  // After a manual action, restart the timer so it doesn't jump right away
  const manual = (action: () => void) => () => {
    action();
    startAutoScroll();
  };

  prevBtn?.addEventListener('click', manual(prev));
  nextBtn?.addEventListener('click', manual(next));
  indicators.forEach((dot, i) => dot.addEventListener('click', manual(() => goTo(i))));

  // Pause while the pointer is over the cards
  carousel.addEventListener('mouseenter', stopAutoScroll);
  carousel.addEventListener('mouseleave', startAutoScroll);

  // Manual scroll (touch, trackpad): follow it with the dots
  let scrollTimeout: ReturnType<typeof setTimeout> | undefined;
  carousel.addEventListener('scroll', () => {
    clearTimeout(scrollTimeout);
    scrollTimeout = setTimeout(() => {
      const index = Math.min(Math.round(carousel.scrollLeft / step()), lastIndex);
      if (index !== currentIndex) {
        currentIndex = index;
        updateIndicators();
      }
    }, 150);
  });

  // Fewer or more cards fit when the window changes size
  let resizeTimeout: ReturnType<typeof setTimeout> | undefined;
  window.addEventListener('resize', () => {
    clearTimeout(resizeTimeout);
    resizeTimeout = setTimeout(() => {
      lastIndex = computeLastIndex();
      goTo(currentIndex, 'auto');
    }, 150);
  });

  lastIndex = computeLastIndex();
  updateIndicators();
  startAutoScroll();
}
