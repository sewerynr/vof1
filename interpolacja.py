import numpy as np

# wywolanie Tinter = inter(wsp_wezl za xy , cells za def_kom, T za T)
# interpolacja na wezly siatki z centrow komorek

def inter(xy, def_kom, T):
    Tn = np.zeros(len(xy))                 # tyle ile punktow wezlow siatki
    d = np.zeros_like(Tn)
    kom_wez = [list() for i in xy]         # tyle ile pkt wezl siatki ale lista list

    for cid, kom in enumerate(def_kom):     # def_kom czyli cells [ 0 1 12 11] itd.
        for w in kom:                       # dla kazdej komorki
            kom_wez[w].append(cid)          # do kom_wez[odp wezlowi siatki] dodaj numer komorki
    for wid, kom in enumerate(kom_wez):
        wxy = xy[wid]                       # pobierz wsp punktu siatki (meszu) w ktorej interpoluje temp
        for c in kom:
            cw = def_kom[c]                  # wezly c1 [0, 1, 12, 11]
            cxy = xy[cw]                     # wsp. wezl dla komorki xy[ [0, 1, 12, 11] ]
            csr = sum(cxy) / len(cxy)        # z nich sr komorki jako suma wsp/ ilosc ws
            delta = wxy - csr                # wsp wektora od sr_kom do
            dist = np.sqrt(np.dot(delta, delta))    # dlugosc tego wektora
            # mam juz wartosc temp w sr kom T[c] i mam odl do punktu siatki w ktorym interp temp
            Tn[wid] += 1./dist * T[c]       # wartosc temp w tym pkt. to
            d[wid] += 1./dist

    Tn = Tn / d
    return Tn


