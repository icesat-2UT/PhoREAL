# -*- mode: python ; coding: utf-8 -*-

import sys
sys.setrecursionlimit(5000)

block_cipher = None


a = Analysis(['getAtl03_GUI.py'],
             pathex=['C:\\Users\\malonzo\\GLAM\\Python\\PhoREAL_GUI\\GUI_3.27'],
             binaries=[],
             datas=[('pho_image.ico', '.')],
             hiddenimports=['pyproj._compat', 'pyproj._datadir', 'pyproj.datadir', 'rasterio.crsrasterio._shim', 'rasterio.compat', 'rasterio.coords', 'rasterio.crs', 'rasterio.dtypes', 'rasterio.enums', 'rasterio.env', 'rasterio.errors', 'rasterio.features', 'rasterio.fill', 'rasterio.mask', 'rasterio.merge', 'rasterio.plot', 'rasterio.profiles', 'rasterio.sample', 'rasterio.tool', 'rasterio.transform', 'rasterio.vfs', 'rasterio.warp', 'rasterio.windows', 'rasterio.__init__'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='PhoReal_v3.27',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True , icon='pho_image.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='PhoReal_v3.27')
