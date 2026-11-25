// Builds the single-file application: index.html and dist/BondX_Studio.html
const fs = require('fs'), path = require('path');
const r = p => fs.readFileSync(path.join(__dirname, p), 'utf8');
const t = r('src/template.html'), ff = r('src/ff_data.js'), eng = r('src/engine.js'), ui = r('src/ui.js'), demo = r('examples/demo_complex.pdb');
for (const [n, s] of [['ff_data', ff], ['engine', eng], ['ui', ui], ['demo', demo]]) if (/<\/script/i.test(s)) throw new Error(n + ' contains </script');
const out = t.replace('__DEMO_PDB__', () => demo).replace('__FFDATA__', () => ff).replace('__ENGINE__', () => eng).replace('__UI__', () => ui);
fs.mkdirSync(path.join(__dirname, 'dist'), { recursive: true });
fs.writeFileSync(path.join(__dirname, 'index.html'), out);
fs.writeFileSync(path.join(__dirname, 'dist', 'BondX_Studio.html'), out);
console.log(`built index.html and dist/BondX_Studio.html (${(out.length / 1024).toFixed(1)} kB)`);
