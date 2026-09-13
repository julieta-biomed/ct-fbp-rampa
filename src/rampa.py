import numpy as np
from skimage.transform import radon, iradon
from fantoma import abdomen, a_mu

def filtro_rampa(sino, tipo='ramp'):
    """Filtra cada proyeccion con |omega| en el dominio de Fourier.

    Esta es la unica diferencia entre retroproyeccion simple y FBP.
    """
    n = sino.shape[0]
    # relleno a potencia de 2 para evitar contaminacion circular
    n_pad = max(64, 2**int(np.ceil(np.log2(2*n))))
    proy = np.zeros((n_pad, sino.shape[1]))
    proy[:n] = sino
    f = np.fft.fftfreq(n_pad).reshape(-1,1)
    rampa = 2*np.abs(f)                       # |omega|, normalizado
    if tipo=='shepp-logan':
        w = np.pi*f/(2*np.abs(f).max()+1e-12)
        rampa = rampa*np.sinc(w/np.pi)
    elif tipo=='hann':
        rampa = rampa*(0.5+0.5*np.cos(np.pi*f/np.abs(f).max()))
    F = np.fft.fft(proy, axis=0)*rampa
    return np.real(np.fft.ifft(F, axis=0))[:n]

def retroproyectar(sino, ang, N):
    """Retroproyeccion pura, sin filtro: 'untar' cada proyeccion de vuelta."""
    return iradon(sino, theta=ang, filter_name=None, circle=True, output_size=N)

if __name__=='__main__':
    N=256; ang=np.linspace(0,180,360,endpoint=False)
    img=abdomen(N, lesion_hu=95); mu=a_mu(img)
    sino=radon(mu, theta=ang, circle=True)

    bp   = retroproyectar(sino, ang, N)
    fbp_manual = retroproyectar(filtro_rampa(sino), ang, N)
    fbp_skimage = iradon(sino, theta=ang, filter_name='ramp', circle=True, output_size=N)

    # escalar a HU para comparar con la verdad
    def a_hu(m, mu_agua=0.02): return 1000*(m/mu_agua - 1)
    mascara = a_mu(img) > 0     # solo dentro del cuerpo
    verdad = img
    print(f'{"metodo":>34} | {"RMSE (HU)":>10} | {"correlacion":>12}')
    print('-'*62)
    for nom, rec in [('retroproyeccion simple', bp),
                     ('FBP con filtro rampa (a mano)', fbp_manual),
                     ('FBP de skimage (ramp)', fbp_skimage)]:
        # ajuste lineal para eliminar escala/offset arbitrarios
        a,b = np.polyfit(rec[mascara].ravel(), verdad[mascara].ravel(), 1)
        est = a*rec+b
        rmse = np.sqrt(((est[mascara]-verdad[mascara])**2).mean())
        cor = np.corrcoef(rec[mascara].ravel(), verdad[mascara].ravel())[0,1]
        print(f'{nom:>34} | {rmse:>10.1f} | {cor:>12.4f}')
    np.save('bp.npy',bp); np.save('fbp.npy',fbp_skimage)
    np.save('fbp_manual.npy',fbp_manual); np.save('img.npy',img); np.save('sino.npy',sino)
