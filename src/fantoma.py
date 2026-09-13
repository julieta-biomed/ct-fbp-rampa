"""Fantoma abdominal con valores de atenuacion en unidades Hounsfield reales."""
import numpy as np

# Valores HU tipicos (Hounsfield 1973; rangos clinicos habituales)
HU = dict(aire=-1000, pulmon=-700, grasa=-90, agua=0, musculo=45,
          higado=60, sangre=55, hueso_esponjoso=300, hueso_cortical=1200,
          metal=3000)

def elipse(N, cy, cx, ry, rx, ang=0.0):
    y,x = np.ogrid[:N,:N]
    y = (y-cy)/ (N/2); x = (x-cx)/(N/2)
    c,s = np.cos(ang), np.sin(ang)
    yr = y*c + x*s; xr = -y*s + x*c
    return (yr/(ry))**2 + (xr/(rx))**2 <= 1

def abdomen(N=256, lesion_hu=None, lesion_r=0.03, metal=False):
    """Corte abdominal simplificado, en unidades Hounsfield.

    lesion_hu: si se indica, inserta una lesion focal en el higado.
    metal: agrega una protesis metalica (para el articulo de artefactos).
    """
    c = N/2
    img = np.full((N,N), float(HU['aire']))
    cuerpo = elipse(N, c, c, 0.86, 0.68)
    img[cuerpo] = HU['grasa']                      # panículo adiposo
    img[elipse(N, c, c, 0.78, 0.60)] = HU['musculo']
    img[elipse(N, c*1.02, c, 0.70, 0.52)] = HU['agua']   # cavidad
    # higado (derecha del paciente = izquierda de la imagen)
    img[elipse(N, c*0.88, c*0.66, 0.30, 0.26, 0.3)] = HU['higado']
    # bazo
    img[elipse(N, c*0.95, c*1.36, 0.16, 0.12, -0.4)] = HU['sangre']
    # vertebra + apofisis
    img[elipse(N, c*1.42, c, 0.13, 0.14)] = HU['hueso_esponjoso']
    img[elipse(N, c*1.42, c, 0.10, 0.11)] = HU['musculo']   # canal medular
    # costillas
    for a in (0.55, 0.95, 1.35, 1.75, 2.25):
        yy = c + 0.72*c*np.cos(a); xx = c + 0.60*c*np.sin(a)
        img[elipse(N, yy, xx, 0.035, 0.055, a)] = HU['hueso_cortical']
        yy2 = c + 0.72*c*np.cos(-a); xx2 = c + 0.60*c*np.sin(-a)
        img[elipse(N, yy2, xx2, 0.035, 0.055, -a)] = HU['hueso_cortical']
    if lesion_hu is not None:
        img[elipse(N, c*0.86, c*0.62, lesion_r, lesion_r)] = lesion_hu
    if metal:
        img[elipse(N, c*1.30, c*0.78, 0.035, 0.035)] = HU['metal']
    return img

def a_mu(hu, mu_agua=0.02):
    """Convierte HU a coeficiente de atenuacion lineal (1/mm).

    HU = 1000 * (mu - mu_agua)/(mu_agua - mu_aire), con mu_aire ~ 0.
    """
    return mu_agua*(1 + hu/1000.0)

def ventana(img, centro, ancho):
    """Aplica una ventana radiologica (window/level) y devuelve 0-1."""
    lo, hi = centro-ancho/2, centro+ancho/2
    return np.clip((img-lo)/(hi-lo), 0, 1)
