import re, csv, glob, sys
H = 0.5e-3
def read_internal(path):
    s = open(path).read()
    k = s.index('internalField')
    m = re.compile(r'nonuniform\s+List<(\w+)>\s*(\d+)\s*\(').search(s, k)
    kind, n, p = m.group(1), int(m.group(2)), m.end()
    if kind == 'scalar':
        return [float(v) for v in s[p:].split(None, n)[:n]]
    out = []
    for mm in re.finditer(r'\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)', s[p:]):
        out.append((float(mm.group(1)), float(mm.group(2)), float(mm.group(3))))
        if len(out) == n: break
    return out
k = lambda v: round(v / H - 0.5)
tdir = sys.argv[1] if len(sys.argv) > 1 else '2'
C = read_internal('0/C'); a = read_internal(tdir + '/alpha.metal')
core = lambda c: -12e-3 < c[0] < 32e-3 and -6e-3 < c[1] < 6e-3
metal = {(k(c[0]), k(c[1]), k(c[2])) for c, av in zip(C, a) if core(c) and av > 0.5}
ev = {}
for f in glob.glob('postProcessing/exacaRecorder/events_p*.csv'):
    for r in list(csv.reader(open(f)))[1:]:
        x, y, z = map(float, r[:3]); key = (k(x), k(y), k(z)); ev[key] = ev.get(key, 0) + 1
x0, x1 = 6e-3, 26e-3   # steady window
bead = {c for c in metal if k(x0) <= c[0] <= k(x1) and c[2] >= 0}
recb = {c for c in ev if k(x0) <= c[0] <= k(x1) and c[2] >= 0}
recs = {c for c in ev if k(x0) <= c[0] <= k(x1) and c[2] < 0}
ns = k(x1) - k(x0) + 1
print(f"steady window x 6..26 mm ({ns} slices of 0.5 mm), metal at t = {tdir} s")
print(f"  bead cells above z=0     : {len(bead):5d}  ({len(bead)/ns:.1f} per slice)")
print(f"  ...with a recorded event : {len(bead & recb):5d}  -> coverage {100*len(bead & recb)/max(len(bead),1):.0f}%")
print(f"  recorded cells below z=0 : {len(recs):5d}  ({len(recs)/ns:.1f} per slice, melted plate)")
kx = k(10.25e-3)
print(f"\n  cross-section at x = 10.25 mm   # metal with event   o metal, no event   . gas   - plate, no event")
for kz in range(6, -6, -1):
    row = ''
    for ky in range(-9, 9):
        c = (kx, ky, kz)
        if c in ev: row += '#'
        elif c in metal: row += 'o' if kz >= 0 else '-'
        else: row += '.'
    print(f"  z {(kz+0.5)*H*1e3:5.2f} mm  {row}" + ("   <- plate top below this row" if kz == 0 else ""))
print(f"  y from {(-9+0.5)*H*1e3:.2f} to {(8+0.5)*H*1e3:.2f} mm, one character = 0.5 mm")
