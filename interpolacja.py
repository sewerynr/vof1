import numpy as np

# wywolanie Tinter = inter(wsp_wezl za xy , cells za def_kom, T za T)

def inter(xy, def_kom, T):
    Tn = np.zeros(len(xy)) # tyle ile wsp_wezl
    d = np.zeros_like(Tn)
    kom_wez = [list() for i in xy]

    for cid, kom in enumerate(def_kom): #def_kom czyli cells [ 0 1 12 11] itd.
        # print cid, kom
        for w in kom:
            kom_wez[w].append(cid)

    for wid, kom in enumerate(kom_wez):
        wxy = xy[wid]
        for c in kom:
            cw = def_kom[c]  # wezly c1 [0, 1, 12, 11]
            cxy = xy[cw] # wsp. wezl c1 xy[ [0, 1, 12, 11] ]
            csr = sum(cxy) / len(cxy)
            delta = wxy - csr
            dist = np.sqrt(np.dot(delta, delta))

            Tn[wid] += 1./dist * T[c]
            d[wid] += 1./dist

    Tn = Tn / d
    return Tn


