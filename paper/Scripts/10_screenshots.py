#!/usr/bin/env python3
"""10_screenshots.py
Interface screenshots for Fig. 3 and the ESI (Playwright + Chromium; the BondX Studio HTML file is opened locally).
The 3D viewer needs the 3Dmol.js library from a CDN; it is hidden when that library cannot be loaded.
Usage: python3 10_screenshots.py /path/to/BondX_Studio.html /path/to/1us0_idd594.pdb ../Figures_V2_MS/shots
"""
import asyncio, os, sys
from playwright.async_api import async_playwright

HTML, PDB, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
HIDE_VIEWER = "document.querySelectorAll('.vslot').forEach(e => { if (!window.$3Dmol) e.style.display = 'none'; });"


async def main():
    os.makedirs(OUT, exist_ok=True)
    async with async_playwright() as p:
        b = await p.chromium.launch(executable_path=os.environ.get('CHROME', '/usr/bin/google-chrome'), args=['--no-sandbox'])
        pg = await b.new_page(viewport={'width': 1440, 'height': 900}, device_scale_factor=2)
        await pg.goto('file://' + os.path.abspath(HTML)); await pg.wait_for_timeout(1200)
        nav = lambda pid: pg.click(f'#side button[data-p={pid}]')
        await pg.set_input_files('#dropPdb input', PDB); await pg.wait_for_timeout(800)
        await pg.evaluate(HIDE_VIEWER); await pg.screenshot(path=os.path.join(OUT, 'a_structure.png'))
        await nav('detect'); await pg.click('#btnDetect'); await pg.wait_for_timeout(500); await pg.evaluate(HIDE_VIEWER)
        await pg.evaluate("document.querySelector('#pg-detect .assume').open = false"); await pg.evaluate('window.scrollTo(0, 0)'); await pg.wait_for_timeout(200)
        await pg.screenshot(path=os.path.join(OUT, 'b_detect.png'))
        await pg.click('#btnToReg'); await nav('qm'); await pg.wait_for_timeout(300)
        await pg.click('#qEnvSeg button[data-env=qmmm]'); await pg.wait_for_timeout(200)
        for bx in await pg.query_selector_all('#qmReg .selc'): await bx.check()
        await pg.wait_for_timeout(500); await pg.evaluate(HIDE_VIEWER)
        await pg.evaluate("document.querySelector('#pg-qm .assume').open = false"); await pg.evaluate('window.scrollTo(0, 0)'); await pg.wait_for_timeout(200)
        await pg.screenshot(path=os.path.join(OUT, 'c_qm_setup.png'))
        await pg.click('#btnBuild'); await pg.wait_for_timeout(1500); await pg.evaluate(HIDE_VIEWER)
        await pg.evaluate("document.querySelector('header.top').style.position = 'static'")
        el = await pg.query_selector('#qOut'); await el.scroll_into_view_if_needed(); await pg.wait_for_timeout(200)
        await el.screenshot(path=os.path.join(OUT, 'd_qm_result.png'))
        await nav('params'); await pg.wait_for_timeout(300); await pg.screenshot(path=os.path.join(OUT, 'e_params.png'))
        await pg.evaluate("document.querySelector('header.top').style.position = 'sticky'")
        await pg.click('#btnTutorial'); await pg.wait_for_timeout(300); await pg.evaluate('window.scrollTo(0, 0)')
        await pg.screenshot(path=os.path.join(OUT, 'f_tutorial_1.png'))
        for i in range(4): await pg.click('#tutNext')
        await pg.wait_for_timeout(300); await pg.evaluate('window.scrollTo(0, 0)'); await pg.screenshot(path=os.path.join(OUT, 'g_tutorial_5.png'))
        await b.close()
    print('screenshots in', OUT)

asyncio.run(main())
