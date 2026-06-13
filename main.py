import pulp
from itertools import combinations, permutations

# =============================================================================
# PARAMETRES
# =============================================================================

annees       = [1, 2, 3, 4, 5]
destinations = ['A', 'C', 'G', 'B', 'H']

dist_liege = {'A': 105, 'C': 100, 'G': 140, 'B': 100, 'H': 60}
dist_entre = {
    ('A', 'C'): 100, ('A', 'G'): 40,  ('A', 'B'): 45,  ('A', 'H'): 50,
    ('C', 'G'): 100, ('C', 'B'): 60,  ('C', 'H'): 80,
    ('G', 'B'): 40,  ('G', 'H'): 60,
    ('B', 'H'): 50
}

def dist(i, j):
    if i == j: return 0
    if (i, j) in dist_entre: return dist_entre[(i, j)]
    if (j, i) in dist_entre: return dist_entre[(j, i)]
    return None

demande_acide = {
    'A': {t: 9000  for t in annees},
    'C': {t: 12000 for t in annees},
    'G': {t: 2000  for t in annees},
    'B': {t: 6200  for t in annees},
    'H': {1: 350, 2: 350, 3: 1300, 4: 1300, 5: 1300},
}
demande_base = 30000

P  = {1: 140000, 2: 200000}
E  = 5000
alpha          = 0.2
cout_carburant = 1.5

CAP_GRAND      = 16.5
CAP_PETIT      =  5.5
H_DISPO        = 1760
VITESSE        = 70
MIN_LIVR       =  5.0
DUREE_RECONFIG = 24

n0 = {1: 4, 2: 6}

dist_base_pure = 2 * dist_liege['A']
tau_base_pure  = dist_base_pure / VITESSE + 1

# =============================================================================
# ENUMERATION DES ROUTES
# =============================================================================

def ordre_optimal(arrets):
    best_dist, best_ordre = float('inf'), list(arrets)
    for perm in permutations(arrets):
        d = dist_liege[perm[0]]
        for i in range(len(perm) - 1):
            d += dist(perm[i], perm[i+1])
        d += dist_liege[perm[-1]]
        if d < best_dist:
            best_dist, best_ordre = d, list(perm)
    return best_ordre, best_dist

routes_acide = {}
for nb in [1, 2, 3]:
    for arrets in combinations(destinations, nb):
        ordre, distance = ordre_optimal(list(arrets))
        tau = distance / VITESSE + nb
        rid = '->'.join(['L'] + ordre + ['L'])
        routes_acide[rid] = {'arrets': ordre, 'distance': distance, 'tau': tau}

routes_mixtes_norm = {}
routes_mixtes_inv  = {}

for rid, rdata in routes_acide.items():
    if 'A' in rdata['arrets']:
        mid = rid + '_MN'
        routes_mixtes_norm[mid] = {
            'arrets_acide': rdata['arrets'],
            'distance':     rdata['distance'],
            'tau':          rdata['tau'],
        }
        mid_inv = rid + '_MI'
        routes_mixtes_inv[mid_inv] = {
            'arrets_acide': rdata['arrets'],
            'distance':     rdata['distance'],
            'tau':          rdata['tau'],
        }
    else:
        arrets_avec_A = rdata['arrets'] + ['A']
        ordre_mix, dist_mix = ordre_optimal(arrets_avec_A)
        tau_mix = dist_mix / VITESSE + len(arrets_avec_A)
        mid = rid + '_via_A_MN'
        routes_mixtes_norm[mid] = {
            'arrets_acide': rdata['arrets'],
            'distance':     dist_mix,
            'tau':          tau_mix,
        }
        mid_inv = rid + '_via_A_MI'
        routes_mixtes_inv[mid_inv] = {
            'arrets_acide': rdata['arrets'],
            'distance':     dist_mix,
            'tau':          tau_mix,
        }

print(f"Routes acide         : {len(routes_acide)}")
print(f"Routes mixtes norm.  : {len(routes_mixtes_norm)}")
print(f"Routes mixtes inv.   : {len(routes_mixtes_inv)}")

# =============================================================================
# MODELE PULP
# =============================================================================

model = pulp.LpProblem("Transport_Chimique_v8", pulp.LpMinimize)

# ---------------------------------------------------------------------------
# Variables de flotte
# ---------------------------------------------------------------------------
n  = {(k,t): pulp.LpVariable(f"n_{k}_{t}",  lowBound=0, cat='Integer') for k in [1,2] for t in annees}
a  = {(k,t): pulp.LpVariable(f"a_{k}_{t}",  lowBound=0, cat='Integer') for k in [1,2] for t in annees}
v  = {(k,t): pulp.LpVariable(f"v_{k}_{t}",  lowBound=0, cat='Integer') for k in [1,2] for t in annees}

n1_acide = {t: pulp.LpVariable(f"n1ac_{t}", lowBound=0, cat='Integer') for t in annees}
n1_base  = {t: pulp.LpVariable(f"n1ba_{t}", lowBound=0, cat='Integer') for t in annees}

n2_ab = {t: pulp.LpVariable(f"n2ab_{t}", lowBound=0, cat='Integer') for t in annees}
n2_ba = {t: pulp.LpVariable(f"n2ba_{t}", lowBound=0, cat='Integer') for t in annees}

delta_1 = {t: pulp.LpVariable(f"delta1_{t}", lowBound=0, cat='Integer') for t in annees}
delta_2 = {t: pulp.LpVariable(f"delta2_{t}", lowBound=0, cat='Integer') for t in annees}
rho1    = {t: pulp.LpVariable(f"rho1_{t}",   lowBound=0, cat='Integer') for t in annees}
rho2    = {t: pulp.LpVariable(f"rho2_{t}",   lowBound=0, cat='Integer') for t in annees}

# ---------------------------------------------------------------------------
# Variables de transport
# ---------------------------------------------------------------------------
T1_r  = {(rid,t): pulp.LpVariable(f"T1_{rid}_{t}", lowBound=0)
          for rid in routes_acide for t in annees}
Q1_rj = {(rid,j,t): pulp.LpVariable(f"Q1_{rid}_{j}_{t}", lowBound=0)
          for rid in routes_acide for j in routes_acide[rid]['arrets'] for t in annees}

T_base_1 = {t: pulp.LpVariable(f"Tbase1_{t}", lowBound=0) for t in annees}
Q_base_1 = {t: pulp.LpVariable(f"Qbase1_{t}", lowBound=0) for t in annees}

TMN_r      = {(mid,t): pulp.LpVariable(f"TMN_{mid}_{t}", lowBound=0)
               for mid in routes_mixtes_norm for t in annees}
QMN_rj     = {(mid,j,t): pulp.LpVariable(f"QMN_{mid}_{j}_{t}", lowBound=0)
               for mid in routes_mixtes_norm for j in routes_mixtes_norm[mid]['arrets_acide'] for t in annees}
QMN_base_r = {(mid,t): pulp.LpVariable(f"QMNbase_{mid}_{t}", lowBound=0)
               for mid in routes_mixtes_norm for t in annees}

TMI_r      = {(mid,t): pulp.LpVariable(f"TMI_{mid}_{t}", lowBound=0)
               for mid in routes_mixtes_inv for t in annees}
QMI_rj     = {(mid,j,t): pulp.LpVariable(f"QMI_{mid}_{j}_{t}", lowBound=0)
               for mid in routes_mixtes_inv for j in routes_mixtes_inv[mid]['arrets_acide'] for t in annees}
QMI_base_r = {(mid,t): pulp.LpVariable(f"QMIbase_{mid}_{t}", lowBound=0)
               for mid in routes_mixtes_inv for t in annees}

# =============================================================================
# FONCTION OBJECTIF
# =============================================================================

cout = []
for t in annees:
    cout.append(P[1]*a[(1,t)] + P[2]*a[(2,t)])
    cout.append(E*(n[(1,t)] + n[(2,t)]))
    for k in [1,2]:
        cout.append(-(P[k] / (1+alpha)**t) * v[(k,t)])
    for rid, rd in routes_acide.items():
        cout.append(cout_carburant * rd['distance'] * T1_r[(rid,t)])
    for mid, md in routes_mixtes_norm.items():
        cout.append(cout_carburant * md['distance'] * TMN_r[(mid,t)])
    for mid, md in routes_mixtes_inv.items():
        cout.append(cout_carburant * md['distance'] * TMI_r[(mid,t)])
    cout.append(cout_carburant * dist_base_pure * T_base_1[t])

model += pulp.lpSum(cout), "Cout_total"

# =============================================================================
# CONTRAINTES
# =============================================================================

# ------------------------------------------------------------------
# C3. Reconfigurations (boucle séparée, AVANT la grande boucle)
# ------------------------------------------------------------------
for t in annees:
    if t == 1:
        model += delta_1[1] == 0, "d1_init"
        model += delta_2[1] == 0, "d2_init"
        model += rho1[1]    == 0, "rho1_init"
        model += rho2[1]    == 0, "rho2_init"
    else:
        model += delta_1[t] >= n1_acide[t] - n1_acide[t-1], f"d1p_{t}"
        model += delta_1[t] >= n1_acide[t-1] - n1_acide[t], f"d1m_{t}"
        model += rho1[t]    == delta_1[t],                   f"rho1_exact_{t}"

        model += delta_2[t] >= n2_ab[t] - n2_ab[t-1],       f"d2p_{t}"
        model += delta_2[t] >= n2_ab[t-1] - n2_ab[t],       f"d2m_{t}"
        model += rho2[t]    == delta_2[t],                   f"rho2_exact_{t}"
# ------------------------------------------------------------------
# Grande boucle : C1, C2, C4, C5, C6, C7, C8
# ------------------------------------------------------------------
for t in annees:

    # C1. Bilan de flotte
    for k in [1,2]:
        n_prev = n0[k] if t == 1 else n[(k,t-1)]
        model += n[(k,t)] == n_prev + a[(k,t)] - v[(k,t)], f"bilan_{k}_{t}"
        model += v[(k,t)] <= n_prev,                        f"vente_max_{k}_{t}"

    # C2. Partition de la flotte
    model += n1_acide[t] + n1_base[t] == n[(1,t)], f"n1split_{t}"
    model += n2_ab[t]    + n2_ba[t]   == n[(2,t)], f"n2split_{t}"
    model += rho1[t] <= n[(1,t)],                   f"rho1_max_{t}"
    model += rho2[t] <= n[(2,t)],                   f"rho2_max_{t}"

    # C4. Satisfaction demande ACIDE
    for j in destinations:
        r1j  = [rid for rid in routes_acide      if j in routes_acide[rid]['arrets']]
        rmnj = [mid for mid in routes_mixtes_norm if j in routes_mixtes_norm[mid]['arrets_acide']]
        rmij = [mid for mid in routes_mixtes_inv  if j in routes_mixtes_inv[mid]['arrets_acide']]
        model += (
            pulp.lpSum(Q1_rj[(rid,j,t)]  for rid in r1j)
          + pulp.lpSum(QMN_rj[(mid,j,t)] for mid in rmnj)
          + pulp.lpSum(QMI_rj[(mid,j,t)] for mid in rmij)
          >= demande_acide[j][t],
            f"demande_acide_{j}_{t}"
        )

    # C5. Satisfaction demande BASE
    model += (
        Q_base_1[t]
      + pulp.lpSum(QMN_base_r[(mid,t)] for mid in routes_mixtes_norm)
      + pulp.lpSum(QMI_base_r[(mid,t)] for mid in routes_mixtes_inv)
      >= demande_base,
        f"demande_base_{t}"
    )

    # C6. Capacités par tournée
    for rid, rd in routes_acide.items():
        model += (pulp.lpSum(Q1_rj[(rid,j,t)] for j in rd['arrets'])
                  <= CAP_GRAND * T1_r[(rid,t)], f"cap_T1_{rid}_{t}")

    for mid, md in routes_mixtes_norm.items():
        model += (pulp.lpSum(QMN_rj[(mid,j,t)] for j in md['arrets_acide'])
                  <= CAP_GRAND * TMN_r[(mid,t)], f"cap_MN_acide_{mid}_{t}")
        model += (QMN_base_r[(mid,t)]
                  <= CAP_PETIT * TMN_r[(mid,t)], f"cap_MN_base_{mid}_{t}")

    for mid, md in routes_mixtes_inv.items():
        model += (QMI_base_r[(mid,t)]
                  <= CAP_GRAND * TMI_r[(mid,t)], f"cap_MI_base_{mid}_{t}")
        model += (pulp.lpSum(QMI_rj[(mid,j,t)] for j in md['arrets_acide'])
                  <= CAP_PETIT * TMI_r[(mid,t)], f"cap_MI_acide_{mid}_{t}")

    model += Q_base_1[t] <= CAP_GRAND * T_base_1[t], f"cap_base1_{t}"

    # C7. Livraison minimale par arrêt
    for rid, rd in routes_acide.items():
        for j in rd['arrets']:
            model += Q1_rj[(rid,j,t)] >= MIN_LIVR * T1_r[(rid,t)], f"minT1_{rid}_{j}_{t}"

    for mid, md in routes_mixtes_norm.items():
        for j in md['arrets_acide']:
            model += QMN_rj[(mid,j,t)] >= MIN_LIVR * TMN_r[(mid,t)], f"minMN_{mid}_{j}_{t}"
        model += QMN_base_r[(mid,t)]   >= MIN_LIVR * TMN_r[(mid,t)], f"min_base_MN_{mid}_{t}"

    for mid, md in routes_mixtes_inv.items():
        for j in md['arrets_acide']:
            model += QMI_rj[(mid,j,t)] >= MIN_LIVR * TMI_r[(mid,t)], f"minMI_{mid}_{j}_{t}"
        model += QMI_base_r[(mid,t)]   >= MIN_LIVR * TMI_r[(mid,t)], f"min_base_MI_{mid}_{t}"

    model += Q_base_1[t] >= MIN_LIVR * T_base_1[t], f"min_base1_{t}"

    # C8. Contraintes de temps
    model += (
        pulp.lpSum(routes_acide[rid]['tau'] * T1_r[(rid,t)] for rid in routes_acide)
        <= H_DISPO * n1_acide[t] - DUREE_RECONFIG * rho1[t],
        f"temps_T1_acide_{t}"
    )
    model += (
        tau_base_pure * T_base_1[t]
        <= H_DISPO * n1_base[t],
        f"temps_T1_base_{t}"
    )
    # Contrainte MN : limitée par les camions en config ab
    model += (
        pulp.lpSum(routes_mixtes_norm[mid]['tau'] * TMN_r[(mid, t)]
                   for mid in routes_mixtes_norm)
        <= H_DISPO * n2_ab[t],
        f"temps_T2_MN_{t}"
    )

    # Contrainte MI : limitée par les camions en config ba
    model += (
        pulp.lpSum(routes_mixtes_inv[mid]['tau'] * TMI_r[(mid, t)]
                   for mid in routes_mixtes_inv)
        <= H_DISPO * n2_ba[t],
        f"temps_T2_MI_{t}"
    )

    # Contrainte globale T2 : somme totale <= capacité totale - pénalité reconfig
    # (conforme à l'équation 21 du rapport, rho2 déduit une seule fois)
    model += (
        pulp.lpSum(routes_mixtes_norm[mid]['tau'] * TMN_r[(mid, t)]
                   for mid in routes_mixtes_norm)
        + pulp.lpSum(routes_mixtes_inv[mid]['tau'] * TMI_r[(mid, t)]
                     for mid in routes_mixtes_inv)
        <= H_DISPO * n[(2, t)] - DUREE_RECONFIG * rho2[t],
        f"temps_T2_global_{t}"
    )

# =============================================================================
# RESOLUTION
# =============================================================================

solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=300)
status = model.solve(solver)

print(f"\n{'='*75}")
print(f"Statut      : {pulp.LpStatus[model.status]}")
if model.status == 1:
    print(f"Coût total  : {pulp.value(model.objective):,.2f} €")
print(f"{'='*75}")

# =============================================================================
# AFFICHAGE DES RESULTATS
# =============================================================================

if model.status != 1:
    print("Pas de solution optimale trouvée.")
else:
    print("\n--- FLOTTE ET AFFECTATION PAR ANNEE ---")
    print(f"{'An':>3} | {'T1':>4} | {'T2':>4} | {'T1-ac':>6} | {'T1-ba':>6} | "
          f"{'T2-ab':>6} | {'T2-ba':>6} | {'rho1':>5} | {'rho2':>5}")
    print("-" * 70)
    for t in annees:
        print(f"{t:>3} | {pulp.value(n[(1,t)]):>4.0f} | {pulp.value(n[(2,t)]):>4.0f} | "
              f"{pulp.value(n1_acide[t]):>6.0f} | {pulp.value(n1_base[t]):>6.0f} | "
              f"{pulp.value(n2_ab[t]):>6.0f} | {pulp.value(n2_ba[t]):>6.0f} | "
              f"{pulp.value(rho1[t]):>5.0f} | {pulp.value(rho2[t]):>5.0f}")

    for t in annees:
        print(f"\n--- TOURNEES ACTIVES (année {t}) ---")

        for rid, rd in routes_acide.items():
            v1 = pulp.value(T1_r[(rid,t)]) or 0
            if v1 > 0.01:
                q = sum(pulp.value(Q1_rj[(rid,j,t)]) or 0 for j in rd['arrets'])
                print(f"  [T1-ac] {rid:35s} | {v1:6.1f} tournées | {q:8.1f}t acide")

        for mid, md in routes_mixtes_norm.items():
            vm = pulp.value(TMN_r[(mid,t)]) or 0
            if vm > 0.01:
                qa = sum(pulp.value(QMN_rj[(mid,j,t)]) or 0 for j in md['arrets_acide'])
                qb = pulp.value(QMN_base_r[(mid,t)]) or 0
                print(f"  [MN]    {mid:35s} | {vm:6.1f} tournées | "
                      f"{qa:8.1f}t acide (grand) + {qb:5.1f}t base (petit)")

        for mid, md in routes_mixtes_inv.items():
            vm = pulp.value(TMI_r[(mid,t)]) or 0
            if vm > 0.01:
                qa = sum(pulp.value(QMI_rj[(mid,j,t)]) or 0 for j in md['arrets_acide'])
                qb = pulp.value(QMI_base_r[(mid,t)]) or 0
                print(f"  [MI]    {mid:35s} | {vm:6.1f} tournées | "
                      f"{qa:8.1f}t acide (petit) + {qb:5.1f}t base (grand)")

        tb1 = pulp.value(T_base_1[t]) or 0
        if tb1 > 0.01:
            print(f"  [B1]    {'L->A(base pure T1)->L':35s} | {tb1:6.1f} tournées | "
                  f"{pulp.value(Q_base_1[t]):8.1f}t base")

    print("\n--- SATISFACTION DES DEMANDES PAR ANNEE ---")
    for t in annees:
        print(f"\n  Année {t} :")
        for j in destinations:
            r1j  = [rid for rid in routes_acide      if j in routes_acide[rid]['arrets']]
            rmnj = [mid for mid in routes_mixtes_norm if j in routes_mixtes_norm[mid]['arrets_acide']]
            rmij = [mid for mid in routes_mixtes_inv  if j in routes_mixtes_inv[mid]['arrets_acide']]
            livre = (
                sum(pulp.value(Q1_rj[(rid,j,t)])  or 0 for rid in r1j)
              + sum(pulp.value(QMN_rj[(mid,j,t)]) or 0 for mid in rmnj)
              + sum(pulp.value(QMI_rj[(mid,j,t)]) or 0 for mid in rmij)
            )
            print(f"    Acide {j}: {livre:8.1f}t / {demande_acide[j][t]}t demandés")

        base_MN   = sum(pulp.value(QMN_base_r[(mid,t)]) or 0 for mid in routes_mixtes_norm)
        base_MI   = sum(pulp.value(QMI_base_r[(mid,t)]) or 0 for mid in routes_mixtes_inv)
        base_pure = pulp.value(Q_base_1[t]) or 0
        total_base = base_MN + base_MI + base_pure
        print(f"    Base    : {total_base:8.1f}t / {demande_base}t demandés")
        print(f"              dont {base_MN:.1f}t mixtes norm. | "
              f"{base_MI:.1f}t mixtes inv. | {base_pure:.1f}t pures T1")

    # Bloc 1 : utilisation des temps
    for t in annees:
        print(f"\nAnnée {t} - utilisation temps :")
        tps_T1_acide = sum(routes_acide[rid]['tau'] * (pulp.value(T1_r[(rid, t)]) or 0)
                           for rid in routes_acide)
        print(f"  T1-acide : {tps_T1_acide:.0f}h / {1760 * int(pulp.value(n1_acide[t]))}h dispo")
        tps_T1_base = tau_base_pure * (pulp.value(T_base_1[t]) or 0)
        print(f"  T1-base  : {tps_T1_base:.0f}h / {1760 * int(pulp.value(n1_base[t]))}h dispo")
        tps_T2_ba = sum(routes_mixtes_inv[mid]['tau'] * (pulp.value(TMI_r[(mid, t)]) or 0)
                        for mid in routes_mixtes_inv)
        print(f"  T2-ba    : {tps_T2_ba:.0f}h / {1760 * int(pulp.value(n2_ba[t]))}h dispo")
        tps_T2_ab = sum(routes_mixtes_norm[mid]['tau'] * (pulp.value(TMN_r[(mid, t)]) or 0)
                        for mid in routes_mixtes_norm)
        print(f"  T2-ab    : {tps_T2_ab:.0f}h / {1760 * int(pulp.value(n2_ab[t]))}h dispo")

    # Bloc 2 : comptage tournées — boucle séparée, même niveau d'indentation que Bloc 1
    print("\n--- COMPTAGE TOURNEES PAR DESTINATION ---")
    for t in annees:
        print(f"\n  Année {t} :")
        for j in destinations:
            total = 0
            for rid, rd in routes_acide.items():
                if j in rd['arrets']:
                    total += pulp.value(T1_r[(rid, t)]) or 0
            for mid, md in routes_mixtes_norm.items():
                if j in md['arrets_acide']:
                    total += pulp.value(TMN_r[(mid, t)]) or 0
            for mid, md in routes_mixtes_inv.items():
                if j in md['arrets_acide']:
                    total += pulp.value(TMI_r[(mid, t)]) or 0
            print(f"    {j} : {total:.0f} tournées")
        total_base = sum(pulp.value(TMI_r[(mid, t)]) or 0 for mid in routes_mixtes_inv)
        total_base += sum(pulp.value(TMN_r[(mid, t)]) or 0 for mid in routes_mixtes_norm)
        total_base += pulp.value(T_base_1[t]) or 0
        print(f"    Base : {total_base:.0f} tournées")