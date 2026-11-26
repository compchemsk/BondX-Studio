// Builds the single-file application: index.html and dist/BondX_Studio.html
const fs = require('fs'), path = require('path');
const r = p => fs.readFileSync(path.join(__dirname, p), 'utf8');
const t = r('src/template.html'), ff = r('src/ff_data.js'), eng = r('src/engine.js'), ui = r('src/ui.js'), demo = r('examples/demo_complex.pdb');
for (const [n, s] of [['ff_data', ff], ['engine', eng], ['ui', ui], ['demo', demo]]) if (/<\/script/i.test(s)) throw new Error(n + ' contains </script');
// Embedded fonts (SIL Open Font License 1.1, see src/fonts/OFL.txt): the page works offline and loads nothing from the web.
const FONTS = [['Instrument Serif', 'InstrumentSerif-Regular', 400, 'normal'], ['Instrument Sans', 'InstrumentSans-Regular', 400, 'normal'],
  ['Instrument Sans', 'InstrumentSans-Italic', 400, 'italic'], ['Instrument Sans', 'InstrumentSans-Bold', 700, 'normal'], ['IBM Plex Mono', 'IBMPlexMono-Regular', 400, 'normal']];
const fontCss = FONTS.map(([fam, file, w, st]) => `@font-face{font-family:"${fam}";src:url(data:font/woff;base64,${fs.readFileSync(path.join(__dirname, 'src/fonts', file + '.woff')).toString('base64')}) format("woff");font-weight:${w};font-style:${st};font-display:swap}`).join('\n');
// Lab logo of the Quantum Simulation of Molecules lab, embedded as a data URI (works offline)
const labLogo = 'data:image/png;base64,' + fs.readFileSync(path.join(__dirname, 'src/qsm_lab_logo.png')).toString('base64');
const out = t.replace('/*__FONTS__*/', () => fontCss).replaceAll('__LABLOGO__', labLogo).replace('__DEMO_PDB__', () => demo).replace('__FFDATA__', () => ff).replace('__ENGINE__', () => eng).replace('__UI__', () => ui);
fs.mkdirSync(path.join(__dirname, 'dist'), { recursive: true });
fs.writeFileSync(path.join(__dirname, 'index.html'), out);
fs.writeFileSync(path.join(__dirname, 'dist', 'BondX_Studio.html'), out);
console.log(`built index.html and dist/BondX_Studio.html (${(out.length / 1024).toFixed(1)} kB)`);
