import numpy as np
import matplotlib.pyplot as plt
import pandas as pd




b1 = float(input('la valeur de b1[m]:'))
db = float(input('la valeur de db[m]:'))
b2 = b1+db
dz = float(input('la valeur de dz[m]:'))
Q = float(input('la valeur de debit Q[m3/s]:'))
ii = float(input('la pente i[%] (EX:3[%]):'))
n = float(input('la valeur de n:'))
kt = float(input('la valeur du coefficient de perte de charge singuliere Kt:'))
q1 = Q / b1
q2 = Q / b2
i = ii/100




if b1 == b2:
        # -----y2^3+(A+dz)*y2^2+B=0
        q = Q / b1
        ycr = ((q ** 2) / 9.81) ** (1 / 3)
        y1 = (q * n / (i ** (0.5))) ** (3 / 5)
        hsmin = ycr + (q ** 2 / (2 * 9.81 * ycr ** 2))
        hs1 = y1 + (q ** 2 / (2 * 9.81 * y1 ** 2))
        hs2 = hs1 - dz

        B = q ** 2 / (2 * 9.81)


        def decrochement():
                global sing, hs2_, remarque, df, y2, y1, hs1, hs2, type_E, celltxt, collumn
                # ----decrochement negative
                if dz < 0:
                        # eq = y2 ** 3 + (-hs1 + dz) * y2 ** 2 + B
                        coef = [1, -hs1 + dz, 0, B]
                        sol = np.roots(coef)
                        if y1 <= ycr:
                                for item in sol:
                                        if (ycr, y1) > item > 0:
                                                y2 = item
                                type_E = "torrentiel"
                        if y1 >= ycr:
                                for item in sol:
                                        if item > y1 and item > ycr:
                                                y2 = item
                                type_E = "fluvial"
                        sing = "decrochement negatif"

                # ----decrochement positive
                hs2_ = hs2
                if dz > 0:
                        if hs2 < hsmin:
                                remarque = "REMARQUE:\n\n-Hs2 est inferieur a Hs(cr)\n-un stockage a été établis (y1++) jusqu'a Hs2=Hs(cr)"
                                hs1 = hsmin + dz
                                hs2 = hsmin
                                # eq = y1 ** 3 + (-hs1) * y1 ** 2 + B
                                coef = [1, -hs1, 0, B]
                                sol = np.roots(coef)
                                l = []
                                for item in sol:
                                        l.append(item)
                                        if 0 < item < max(l):
                                                y1 = item
                                        else:
                                                y1 = max(l)

                        # eq = y2 ** 3 + (-hs1 - dz) * y2 ** 2 + B
                        coef = [1, -hs1 + dz, 0, B]
                        sol = np.roots(coef)

                        if y1 <= ycr:
                                for item in sol:
                                        if ycr > item > y1:
                                                y2 = item
                                        else:
                                                y2 = ycr
                                type_E = "torrentiel"

                        if y1 >= ycr:
                                for item in sol:
                                        if y1 > item > ycr:
                                                y2 = item
                                        else:
                                                y2 = ycr
                                type_E = "fluvial"
                        sing = "decrochement positif"
                if dz == 0:
                        y2 = y1
                        type_E = "stable"
                        sing = "Null"
                # -----tableau pour les commentaires---
                data = {"data": ["type d'écoulement", "Q[m3/s]", "B1[m]", "B2[m]", "Delta(Z)[m]", "type de singularité",
                                 "i", "n", "Hs(cr)[m]", "Hs(1)[m]", "Hs(2)[m]", "y(cr)[m]", "y(1)[m]", "y(2)[m]"],
                        "resultat": [type_E, Q, b1, b2, dz, sing, i, n, hsmin, hs1, hs2, ycr, y1, y2]}
                df = pd.DataFrame(data)


        def tracage_d():
                decrochement()

                def equation(y):
                        return y + B / (y ** 2)

                # Générer des valeurs y
                y = np.linspace(0.00001, 10, 10000)

                # Calculer les valeurs hs pour chaque y
                hs = equation(y)
                d = y
                # Tracer la courbe
                fig, (f, t) = plt.subplots(1, 2)
                plt.plot(hs, y, label='y = f(Hs)')
                plt.plot(hs1, y1, label='point 1', marker='o', color='black')
                plt.plot(hsmin, ycr, label='point critique', marker='o', color='red')
                if hs2 == hsmin:
                        plt.plot(hs2, y2, label='point 2', marker='o', color='red')
                else:
                        plt.plot(hs2, y2, label='point 2', marker='o', color='blue')

                # tracage du tableau d'information

                table = pd.plotting.table(f, df, loc="center left")
                table.auto_set_font_size(False)
                table.set_fontsize(7)
                f.axis('off')
                plt.subplots_adjust(wspace=0.3)
                if hs2 != hs2_:
                        plt.figtext(x=0.12, y=0.1, s=remarque).set_fontsize(10)

                plt.plot(d, y)
                plt.axhline(0, color='black', linewidth=0.5)
                plt.axvline(0, color='black', linewidth=0.5)
                plt.ylim(0, 10)
                plt.xlim(0, 10)
                plt.title('y = f(Hs)')
                plt.xlabel('Hs')
                plt.ylabel('y')
                plt.legend()
                plt.grid(color='gray', linestyle='--', linewidth=0.5)
                plt.show()
        tracage_d()
#-------ELARGISSEMENT +-DZ--------
if b2 > b1:
        y1 = (q1 * n / (i ** (0.5))) ** (3 / 5)
        hs1b = y1 + ((q1 ** 2) / (2 * 9.81 * y1 ** 2)) * (1 + kt)
        B = ((q2 ** 2) / (2 * 9.81)) * (1 + kt)


        def elargissement():
                global sing, dz, hs1, hs2, h1s, h2s, type_E, df, dz, y2, y1, hs1b, hs2b, hsb2, hsbcr, ycr, d, hsb1, y
                # -----type de singularité-----
                sing = 'elargissement'
                if 0.1 < kt < 0.2:
                        sing = 'elar progressif'
                if kt == 1:
                        sing = 'elar vifs'
                if dz > 0:
                        sing = sing + ' + ' + 'dec (+)'
                if dz < 0:
                        sing = sing + ' + ' + 'dec (-)'

                # ----- eq = y2 ** 3 + (-hs1 - kt*q1**2/2gy1**2) * y2 ** 2 + B
                coef = [1, -dz - hs1b, 0, B]
                s = np.roots(coef)
                sol = np.real(s)

                def equation1(y):
                        return y + (kt + 1) * (q1 ** 2) / (2 * 9.81 * y ** 2)

                def equation2(y):
                        return y + (kt + 1) * (q2 ** 2) / (2 * 9.81 * y ** 2)

                def equation1_(y):
                        return y + (q1 ** 2) / (2 * 9.81 * y ** 2)

                def equation2_(y):
                        return y + (q2 ** 2) / (2 * 9.81 * y ** 2)

                # ----- Générer des valeurs y

                y = np.linspace(0.00001, 10, 10000)

                # ----- Calculer les valeurs hs pour chaque y
                hsb1 = equation1(y)
                hsb2 = equation2(y)
                h1s = equation1_(y)
                h2s = equation2_(y)
                r = 0
                for item in hsb1:
                        if item == equation1(y).min():
                                hsbcr = item
                                ycr = y[r]
                        r += 1
                d = y
                if y1 > ycr > 0:
                        for item in sol:
                                if item > y1:
                                        y2 = item
                        type_E = "fluvial"
                if ycr > y1 > 0:
                        for item in sol:
                                if y1 > item > 0:
                                        y2 = item
                        type_E = "torrentiel"
                hs2b = y2 + ((q2 ** 2) / (2 * 9.81 * y2 ** 2)) * (1 + kt)
                hs1 = y1 + ((q1 ** 2) / (2 * 9.81 * y1 ** 2))
                hs2 = y2 + ((q2 ** 2) / (2 * 9.81 * y2 ** 2))
                dh = (hs1b - hs1) - (hs2b - hs2)
                data = {"data": ["type d'écoulement", "Q[m3/s]", "B1[m]", "B2[m]", "Delta(Z)[m]", "i", "n", "Kt",
                                 "type du singularité", "Hsb(cr)[m]",
                                 "Hsb(1)[m]", "Hsb(2)[m]",
                                 "Hs1[m]", "Hs2[m]", "dh", "y(cr)[m]", "y(1)[m]", "y(2)[m]"],
                        "resultat": [type_E, Q, b1, b2, dz, i, n, kt, sing, equation2(y).min(), hs1b, hs2b, hs1, hs2,
                                     dh, ycr,
                                     y1, y2]}
                df = pd.DataFrame(data)


        def tracage_e():
                elargissement()
                # Tracage de la courbe
                fig, (f, t, h) = plt.subplots(1, 3, sharex=True)

                # la fonction hs1
                t.plot(h1s, y, label='q1', color="green")
                t.set_xlabel('Hs')
                t.set_ylabel('y')
                t.set_title('y = f(hs)')
                t.grid(color='gray', linestyle='--', linewidth=0.5)
                # la fonction hsb1
                h.plot(hsb1, y, label='q1', color="green")
                # les coordonnees du 1er point
                h.plot(hs1b, y1, label='y1,hsb1', marker='o', color="green")
                t.plot(hs1, y1, label='y1,hs1', marker='o', color="green")
                t.legend()
                # la fonction hs2
                t.plot(h2s, y, label='q2', color="blue")
                # la fonction hsb2
                h.plot(hsb2, y, label='q2', color="blue")
                # les coordonnees du 2eme point
                h.plot(hs2b, y2, label='y2,hsb2', marker='o', color="blue")
                t.plot(hs2, y2, label='y2,hs2', marker='o', color="blue")
                t.legend()
                # la fonction y = hs-bar-
                h.plot(d, y)
                # la fonction y = hs
                t.plot(d, y, color="orange")
                # les coordonnees du point critique
                h.plot(hsbcr, ycr, label='point critique', marker='o', color="red")

                # tracage du tableau d'information

                table = pd.plotting.table(f, df, loc="center left")
                table.auto_set_font_size(False)
                table.set_fontsize(7)
                f.axis('off')
                plt.subplots_adjust(wspace=0.3)

                # grillage du graph
                plt.axhline(0, color='black', linewidth=0.5)
                plt.axvline(0, color='black', linewidth=0.5)
                plt.ylim(0, 10)
                plt.xlim(0, 10)
                plt.title("y = f(Hs-bar-)")
                plt.xlabel('Hs-bar-')
                plt.ylabel('y')
                plt.legend()
                plt.grid(color='gray', linestyle='--', linewidth=0.5)
                plt.show()
        tracage_e()
#-------RETRECISSEMENT +-DZ--------
if b1 > b2:
        y2 = (q2 * n / (i ** (0.5))) ** (3 / 5)
        hs2b = y2 + ((q2 ** 2) / (2 * 9.81 * y2 ** 2)) * (1 + kt)
        B = ((q1 ** 2) / (2 * 9.81)) * (1 + kt)


        def retrecissement():
                global  sing, dz, hs1, hs2, h1s, h2s, type_E, df, dz, y2, y1, hs1b, hs2b, hsb2, hsbcr, ycr, d, hsb1, y
                sing = 'retrecissement'
                if 0.05 < kt < 0.1 and db < 0:
                        sing = 'ret progressif'
                if kt == 0.5 and db < 0:
                        sing = 'ret vifs'
                if dz > 0:
                        sing = sing + ' + ' + 'dec (+)'
                if dz < 0:
                        sing = sing + ' + ' + 'dec (-)'

                # eq = y1 ** 3 + (-hs2 - kt*q2**2/2gy2**2) * y1 ** 2 + B
                coef = [1, -dz - hs2b, 0, B]
                s = np.roots(coef)
                sol = np.real(s)

                def equation1(y):
                        return y + (kt + 1) * (q1 ** 2) / (2 * 9.81 * y ** 2)

                def equation2(y):
                        return y + (kt + 1) * (q2 ** 2) / (2 * 9.81 * y ** 2)

                def equation1_(y):
                        return y + (q1 ** 2) / (2 * 9.81 * y ** 2)

                def equation2_(y):
                        return y + (q2 ** 2) / (2 * 9.81 * y ** 2)

                # Générer des valeurs y

                y = np.linspace(0.00001, 10, 10000)

                # Calculer les valeurs hs pour chaque y
                hsb1 = equation1(y)
                hsb2 = equation2(y)
                h1s = equation1_(y)
                h2s = equation2_(y)
                r = 0
                for item in hsb2:
                        if item == equation2(y).min():
                                hsbcr = item
                                ycr = y[r]
                        r += 1
                d = y
                if y2 > ycr > 0:
                        for item in sol:
                                if y2 < item:
                                        y1 = item
                        type_E = "fluvial"
                if ycr > y2 > 0:
                        for item in sol:
                                if y2 > item > 0:
                                        y1 = item
                        type_E = "torrentiel"
                bcr = ((Q ** 2) / (9.81 * (ycr ** 3))) ** 0.5
                if b2 == bcr:
                        y2 = ycr
                        hs2b = hsbcr
                hs1b = y1 + ((q1 ** 2) / (2 * 9.81 * y1 ** 2)) * (1 + kt)
                hs1 = y1 + ((q1 ** 2) / (2 * 9.81 * y1 ** 2))
                hs2 = y2 + ((q2 ** 2) / (2 * 9.81 * y2 ** 2))
                dh = (hs1b - hs1) - (hs2b - hs2)
                data = {"data": ["type d'écoulement", "Q[m3/s]", "B1[m]", "B2[m]", "Delta(Z)[m]", "i", "n", "Kt",
                                 "type du singularité", "Hsb(cr)[m]",
                                 "Hsb(1)[m]", "Hsb(2)[m]",
                                 "Hs1[m]", "Hs2[m]", "dh", "y(cr)[m]", "y(1)[m]", "y(2)[m]"],
                        "resultat": [type_E, Q, b1, b2, dz, i, n, kt, sing, equation2(y).min(), hs1b, hs2b, hs1, hs2,
                                     dh, ycr,
                                     y1, y2]}
                df = pd.DataFrame(data)


        def tracage_r():
                retrecissement()
                # Tracer la courbe
                fig, (f, t, h) = plt.subplots(1, 3, sharex=True)

                # la fonction hs1
                t.plot(h1s, y, label='q1', color="green")
                t.set_xlabel('Hs')
                t.set_ylabel('y')
                t.set_title('y = f(hs)')
                t.grid(color='gray', linestyle='--', linewidth=0.5)

                # la fonction hsb1
                h.plot(hsb1, y, label='q1', color="green")
                # les coordonnees du 1er point
                h.plot(hs1b, y1, label='y1,hsb1', marker='o', color="green")
                t.plot(hs1, y1, label='y1,hs1', marker='o', color="green")
                t.legend()
                # la fonction hs2
                t.plot(h2s, y, label='q2', color="blue")
                # la fonction hsb2
                h.plot(hsb2, y, label='q2', color="blue")
                # les coordonnees du 2eme point
                h.plot(hs2b, y2, label='y2,hsb2', marker='o', color="blue")
                t.plot(hs2, y2, label='y2,hs2', marker='o', color="blue")
                t.legend()
                # la fonction y = hs-bar-
                h.plot(d, y, color="orange")
                # la fonction y = hs
                t.plot(d, y, color="orange")
                # les coordonnees du point critique
                h.plot(hsbcr, ycr, label='point critique', marker='o', color="red")
                # tracage du tableau d'information

                table = pd.plotting.table(f, df, loc="center left")
                table.auto_set_font_size(False)
                table.set_fontsize(7)
                f.axis('off')
                plt.subplots_adjust(wspace=0.3)

                # grillage du graph
                plt.axhline(0, color='black', linewidth=0.5)
                plt.axvline(0, color='black', linewidth=0.5)
                plt.ylim(0, 10)
                plt.xlim(0, 10)
                plt.title("y = f(Hs-bar-)")
                plt.xlabel('Hs-bar-')
                plt.ylabel('y')
                plt.legend()
                plt.grid(color='gray', linestyle='--', linewidth=0.5)
                plt.show()
        tracage_r()
