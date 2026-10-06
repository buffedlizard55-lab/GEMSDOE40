'use strict';
document.querySelectorAll('[data-copy]').forEach(button => {
  button.addEventListener('click', async () => {
    const value = document.getElementById(button.dataset.copy)?.textContent || '';
    try {
      await navigator.clipboard.writeText(value.trim());
      button.textContent = 'Copied';
    } catch (_) {
      button.textContent = 'Select the text above to copy';
    }
  });
});
const search = document.getElementById('project-search');
if (search) search.addEventListener('input', () => {
  const query = search.value.toLowerCase();
  let visible = 0;
  document.querySelectorAll('[data-project-row]').forEach(row => {
    const show = row.textContent.toLowerCase().includes(query);
    row.classList.toggle('hide', !show);
    if (show) visible += 1;
  });
  const count = document.getElementById('search-count');
  if (count) count.textContent = `${visible} project rows shown`;
});
async function contextFeed() {
  const status = document.getElementById('feed-status');
  if (!status) return;
  try {
    const response = await fetch('data/source-feed.json', {cache: 'no-cache'});
    if (!response.ok) throw new Error('feed unavailable');
    const feed = await response.json();
    const time = feed.last_success_utc ? new Date(feed.last_success_utc).toLocaleString(undefined, {timeZone: 'UTC'}) + ' UTC' : 'no successful refresh yet';
    const elapsed = feed.last_success_utc ? Date.now() - new Date(feed.last_success_utc).getTime() : Infinity;
    const stale = feed.status !== 'ok' || elapsed > 48 * 3600 * 1000;
    status.textContent = `${stale ? 'STALE / UNAVAILABLE' : 'SOURCE SNAPSHOT'} · ${time}. Daily build refresh; not a real-time fault detector.${feed.source_probe_failures ? ` ${feed.source_probe_failures} metadata probe(s) unavailable; see source JSON.` : ''}`;
    const list = document.getElementById('feed-events');
    if (!list) return;
    list.replaceChildren();
    (feed.events || []).slice(0, 5).forEach(event => {
      const li = document.createElement('li');
      const a = document.createElement('a');
      const url = new URL(event.url);
      if (url.protocol !== 'https:' || url.hostname !== 'earthquake.usgs.gov') return;
      a.href = url.href;
      a.textContent = `M ${event.magnitude ?? '?'} · ${event.place || 'USGS event'}`;
      li.appendChild(a);
      list.appendChild(li);
    });
    if (!list.children.length) {
      const li = document.createElement('li');
      li.textContent = feed.status === 'ok' ? 'USGS returned no matching events for this interval.' : 'No event snapshot available. Open the official USGS source below; no data are fabricated.';
      list.appendChild(li);
    }
  } catch (_) {
    status.textContent = 'Feed unavailable. The verified GeoTIFF remains downloadable; source status does not change the submission gate.';
  }
}
contextFeed();
