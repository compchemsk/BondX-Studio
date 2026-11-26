#!/usr/bin/env python3
"""23_v3_screenshots.py - interface screenshots of BondX Studio 3.0 for Fig. 2 and Fig. S11 (Playwright + Chromium).
Usage: python3 23_v3_screenshots.py <BondX_Studio.html> <paper folder>
For the figures all page text is set in Arial (Liberation Sans where Arial is not installed); the app itself uses its own fonts.
The 3D viewer needs 3Dmol.js from a CDN; offline it is hidden. Panel (c) reads a real NBO output of the paper (1US0 cluster model
of V2, with its manifest), so it shows the real decision values."""
import asyncio, os, sys
from playwright.async_api import async_playwright
HTML, ROOT = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
OUT = os.path.join(ROOT, 'MS_Figures_V3', 'shots'); os.makedirs(OUT, exist_ok=True)
PDB = os.path.join(ROOT, 'data', 'protein', '1us0_idd594.pdb')
MAN = os.path.join(ROOT, 'qm_package', 'F_protein', 'manifests', '1us0_idd594__LDTA__cluster_pcm4_manifest.json')
LOG = os.path.join(ROOT, 'qm_package', 'F_protein', 'gaussian', '1us0_idd594__LDTA__cluster_pcm4.log')
ARIAL = ".toast{display:none !important} *{font-family:Arial,'Liberation Sans',sans-serif !important} .pagehead h2,.brand h1,.card h2,.tuth{font-family:Arial,'Liberation Sans',sans-serif !important;font-weight:700 !important;letter-spacing:0 !important} .mono,code,pre{font-family:'Liberation Mono',monospace !important}"
HIDE = "document.querySelectorAll('.vslot').forEach(e => { if (!window.$3Dmol) e.style.display = 'none'; }); document.querySelectorAll('.assume').forEach(e => e.open = false);"

async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(); pg = await b.new_page(viewport={'width': 1600, 'height': 960}, device_scale_factor=2)
        await pg.goto('file://' + HTML); await pg.wait_for_timeout(1500)
        await pg.evaluate("localStorage.setItem('bondx-theme','light')"); await pg.reload(); await pg.wait_for_timeout(1500)
        await pg.add_style_tag(content=ARIAL + " .reveal-on .page.on .pagehead,.reveal-on .page.on .card,.reveal-on .page.on .assume{opacity:1 !important;transform:none !important}")
        nav = lambda pid: pg.click(f'#side button[data-p={pid}]')
        async def shot(name, sel=None):
            await pg.evaluate(HIDE); await pg.evaluate('window.scrollTo(0,0)'); await pg.wait_for_timeout(400)
            if sel: el = await pg.query_selector(sel); await el.screenshot(path=os.path.join(OUT, name + '.png'))
            else: await pg.screenshot(path=os.path.join(OUT, name + '.png'))
        async def clip(name, page, first, n):
            await pg.evaluate(HIDE); await pg.wait_for_timeout(300)
            bb = await pg.evaluate(f"""(() => {{ const kids = [...document.querySelector('{page}').children].filter(e => !e.classList.contains('pagehead') && !e.classList.contains('assume') && e.offsetHeight > 0);
              const sel = kids.slice({first}, {first} + {n}); const r = sel.map(e => e.getBoundingClientRect()); const top = Math.min(...r.map(x => x.top)) + scrollY, bot = Math.max(...r.map(x => x.bottom)) + scrollY;
              const left = Math.min(...r.map(x => x.left)), right = Math.max(...r.map(x => x.right)); return [left - 8, top - 8, right - left + 16, bot - top + 16]; }})()""")
            await pg.screenshot(path=os.path.join(OUT, name + '.png'), full_page=True, clip=dict(x=bb[0], y=bb[1], width=bb[2], height=bb[3]))
        await pg.set_input_files('#dropPdb input', PDB); await pg.wait_for_timeout(900); await shot('s01_structure')
        await nav('detect'); await pg.click('#btnDetect'); await pg.wait_for_timeout(600); await shot('a_detect'); await shot('s02_detect'); await clip('a_detect_clip', '#pg-detect', 0, 2)
        await pg.click('#btnToReg'); await nav('qm'); await pg.wait_for_timeout(300)
        await pg.click('#qEnvSeg button[data-env=qmmm]'); await pg.wait_for_timeout(200)
        for bx in await pg.query_selector_all('#qmReg .selc'): await bx.check()
        await pg.wait_for_timeout(600); await shot('b_qm_setup'); await clip('b_qm_clip', '#pg-qm', 0, 2)
        await pg.click('#btnBuild'); await pg.wait_for_timeout(1800)
        await pg.evaluate("document.querySelector('header.top').style.position='static'")
        el = await pg.query_selector('#qOut'); await el.scroll_into_view_if_needed(); await pg.wait_for_timeout(300); await el.screenshot(path=os.path.join(OUT, 's03_qm_result.png'))
        await pg.evaluate("document.querySelector('header.top').style.position=''")
        await nav('nbo'); await pg.wait_for_timeout(300)
        await pg.set_input_files('#dropMan input', MAN); await pg.wait_for_timeout(500)
        await pg.set_input_files('#dropNbo input', LOG); await pg.wait_for_timeout(1200); await shot('c_nbo'); await shot('s04_nbo')
        card = pg.locator('#nRes').locator('xpath=ancestor::div[contains(@class,"card")][1]'); await card.screenshot(path=os.path.join(OUT, 'c_nbo_card.png'))
        await nav('params'); await pg.wait_for_timeout(300); await shot('s05_params')
        await nav('method'); await pg.wait_for_timeout(300); await shot('s06_method')
        await pg.click('#btnTutorial'); await pg.wait_for_timeout(300); await shot('s07_tutorial1')
        for i in range(4): await pg.click('#tutNext')
        await pg.wait_for_timeout(300); await shot('s08_tutorial5')
        await b.close()
    print('screenshots in', OUT)
asyncio.run(main())
