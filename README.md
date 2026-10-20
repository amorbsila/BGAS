**Zenodo DOI:** [10.5281/zenodo.23518746](https://doi.org/10.5281/zenodo.23518746)

**Archived release:** [v1.0.1](https://github.com/amorbsila/BGAS/releases/tag/v1.0.1)
# BGAS — Reconstruction numérique en base de Pell

**Auteur : Amor Bsila**

## 1. Présentation du projet

Ce projet étudie le système BGAS de reconstruction numérique dans la base

\[
\rho=1+\sqrt{2},
\qquad
\bar{\rho}=1-\sqrt{2}.
\]

L'objectif est d'examiner les relations entre les codes binaires BGAS, les nombres de Pell, les nombres de Pell–Lucas et certaines structures arithmétiques associées.

Le travail associe des formules mathématiques, des algorithmes Python et des vérifications numériques reproductibles.

## 2. Principe de reconstruction BGAS

Pour un mot binaire \(w=(w_1,\ldots,w_k)\), la reconstruction `BGASINV` part de \(u_0=1\) et applique successivement

\[
u_j=\operatorname{round}(\rho u_{j-1})+w_j,
\]

où l'arrondi est effectué à l'entier le plus proche, avec la convention *ROUND_HALF_UP*.

L'étude porte notamment sur les mots constitués de zéros, de uns et sur des motifs binaires particuliers.

## 3. Relations avec les suites de Pell

Les nombres de Pell sont définis par

\[
P_0=0,\quad P_1=1,\quad
P_n=2P_{n-1}+P_{n-2}.
\]

Leurs premiers termes sont

\[
0,1,2,5,12,29,70,169,\ldots
\]

Les nombres de Pell–Lucas sont étudiés à travers les deux suites

\[
H_n=\frac{\rho^n+\bar{\rho}^{\,n}}{2},
\qquad Q_n=2H_n.
\]

Les relations BGAS examinées comprennent notamment :

\[
\operatorname{BGASINV}(0^\ell)=P_{\ell+1},
\]

\[
\operatorname{BGASINV}(10^{\ell-2})=H_\ell,
\]

\[
\operatorname{BGASINV}(010^{\ell-2})=Q_\ell.
\]

Ces relations sont étudiées avec leurs conventions d'indices et les conditions indiquées dans la note mathématique.

## 4. Contenu du dépôt

- **`bgas_outils.py`** — outils Python utilisés pour les calculs BGAS.
- **`verif_pell.py`** — programme de vérification des formules et des relations étudiées.
- **`note_pell.tex`** — note mathématique présentant les définitions, formules, motifs et observations de l'étude.

## 5. Exécuter les vérifications

Placer les fichiers Python dans le même répertoire, puis exécuter :

```bash
python3 verif_pell.py
```

Le script importe `bgas_outils.py` et effectue les tests prévus dans le programme.

Les résultats dépendent des tests effectivement présents dans la version du script utilisée. Une vérification numérique sur un ensemble fini d'indices ne constitue pas, à elle seule, une démonstration générale.

## 6. Démarche scientifique

Le projet distingue trois niveaux :

1. **Identités mathématiques** : relations établies par une démonstration.
2. **Vérifications numériques** : égalités contrôlées par le programme sur des plages d'indices données.
3. **Conjectures et observations structurelles** : phénomènes observés qui demandent encore une preuve générale.

Cette distinction est particulièrement importante dans l'étude des motifs BGAS, des fibres de codes, des intervalles de Pell et des structures associées aux nombres de Pell–Lucas.

## 7. Reproductibilité et évolution

Les scripts permettent de reproduire les calculs à partir des fichiers publiés. Les modifications et les nouveaux résultats doivent être accompagnés, autant que possible, des tests correspondants.

Le dépôt a vocation à conserver une trace du développement de l'étude : définitions, résultats, vérifications et questions encore ouvertes.

## 8. Licence

Aucune licence de réutilisation n'est déclarée dans cette version. Les conditions de réutilisation pourront être précisées séparément par l'auteur.
