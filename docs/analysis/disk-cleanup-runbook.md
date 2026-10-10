# Disk Cleanup Runbook — WSL case tree and the Windows host

Written 2026-09-18, after the host drive reached 3 GB free while M3 was running.
Updated 2026-09-19: the WSL distro was moved from C: to D: (§1), and the
`--set-sparse` advice was withdrawn (§4, Phase 3).
Reusable: this is the procedure for every future cleanup, not a one-off.

---

## 1. The machine, measured

| | value |
|---|---|
| physical disk | **WD_BLACK SN7100 500 GB NVMe** (Disk 0), split into C: and D: |
| Windows drive C: | 200 GB partition. 3.0 GB free on 2026-09-18 → ~17 GB after uninstalls → ~68 GB after the WSL move |
| drive D: | 265 GB partition, same NVMe. 70 GB free before the move → ~19 GB after |
| WSL distro | Ubuntu-22.04, WSL 2.7.14.0, kernel 6.18.33.2-2 |
| **ext4.vhdx path** | **`D:\WSL\Ubuntu-22.04\ext4.vhdx`** (since 2026-09-19) |
| old path (deleted) | `C:\Users\User\AppData\Local\wsl\{ab3f8a02-dcfe-4027-a361-272931ccf321}\ext4.vhdx` |
| ext4.vhdx size | 51.2 GB at the move (46.6 GB on 2026-09-18) |
| **sparse flag** | **NOT set, and do not set it** (Phase 3) |
| inside WSL | 1007 GB virtual. The 1007 is meaningless: it is the vhdx *maximum* size |
| pagefile.sys | 19.8 GB at C:\ root |
| `~/thesis` | 25 GB, of which `cases` is 25 GB |
| `~/ISW_CFD_CA` | 15 GB, not in the thesis repo |

### The move to D: (2026-09-19), for the record

- D: was confirmed as a partition of the same internal NVMe, so there is no
  speed penalty. `Get-Partition` showed C and D both on DiskNumber 0, and
  `Get-Disk` showed BusType NVMe.
- There was no `wsl --export` backup, because there was no room: the tar
  (~45 GB) plus the moved vhdx (51.2 GB) exceeded the 70 GB free on D:. Code,
  case recipes and results were already on GitHub (tag `0p5mm-M1-M4`), and
  `--move` copies before it deletes.
- Commands, from a Windows PowerShell:
  `wsl --shutdown` → `wsl --list --verbose` (Stopped) →
  `wsl --manage Ubuntu-22.04 --move D:\WSL\Ubuntu-22.04`, which printed
  "The operation completed successfully."
- Verified: registry `BasePath` = `D:\WSL\Ubuntu-22.04`, and the old C: file
  is gone (`Test-Path` = False). `--move` keeps the registry entry, so the
  default user is unchanged, unlike export/import.

### The trap that governs everything here

`df -h` **inside WSL is not the drive.** It reports the vhdx's 1007 GB maximum
size, which has nothing to do with free space on the Windows drive that holds
the vhdx (D: since the move). In the 2026-09-18 incident it showed 911 GB free
while Windows had 3 GB.

Worse, the vhdx is **not sparse**. It grows monotonically and **never shrinks
on its own**. Deleting files inside WSL returns space to the ext4 filesystem
but returns *nothing* to Windows until the vhdx is explicitly compacted.

So deleting inside WSL has two separate effects, and it is worth being clear
about which one you are buying:

1. **Immediate, no restart:** freed ext4 blocks are reused by subsequent
   writes instead of the vhdx expanding. This is what protects a running
   simulation when the host drive is nearly full. Real, but invisible to
   Windows.
2. **Deferred, needs the distro stopped:** compaction actually hands the bytes
   back to Windows.

With ~19 GB free on D: after the move, (1) is the headroom that matters.
Production runs move to the remote PC, which should keep growth here modest.

---

## 2. What a case directory is made of

An OpenFOAM case holds two things with opposite value-per-byte.

**The recipe:** `system/`, `constant/` (minus `polyMesh`), `0.orig/`,
`Allrun`, `Allclean`, `caseParameters`. A few hundred kB. This *is* the run:
mesh spec, schemes, solver tolerances, physical properties, boundary
conditions, arc and droplet parameters. `blockMesh` + `setFields` + the solver
regenerate everything else from it deterministically. Since 2026-09-19 the
recipes are on GitHub.

**The fields:** every numeric time directory, `constant/polyMesh`,
`processor*/`, `log*`. Gigabytes. Expensive in CPU hours, but a deterministic
function of the recipe and the solver binary. They are **not** on GitHub.

The project's standing rule is *"a reference run's fields are the asset;
summary text is a derivative. Text regenerates from fields in seconds; fields
never regenerate from text."* That rule still holds for reference runs. Its
corollary is what makes this cleanup safe: **for a superseded run, the recipe
is the asset and the fields are the derivative**, because nothing will be
measured from those fields again. Archive the recipe, drop the fields.

---

## 3. Classification

### KEEP — never slimmed, never purged

| case | why |
|---|---|
| `waam-arc-01-P1114W` | gate-ON reference; every M-series delta is measured against it |
| `waam-arc-01-M1-mushy` | M1 result: liquidus gate removed (+33.6% pool U, width null) |
| `waam-arc-01-M2-full` | M2 result: Lorentz + droplet-restart + darcy_q (w +18.1%). The 0.5 mm level of the mesh study, and the ExaCA pipeline test case |
| `waam-arc-01-M3-darcyA1e5` | M3: darcy_A 1e5, the creep artifact |
| `waam-arc-01-M4-1512W` | M4: literal η = 0.70, melts into the plate underside |
| `hartmann1d-postpatch` | **Gate A**: u(0) 0.0158753, Ha 15.810, err 1.3869e-02 |
| `hartmann1d` | Gate A unpatched control, 3 MB, proves the patch is a no-op |
| `Plate2D-mhd-postpatch` | **Gate B**: 1869.0 / 0.0504 / 3.4135 / 0.0552 / 12.7787 |
| `Plate2D-mhd` | Gate B reference case |

The M-runs are thesis results in their own right. A later run superseding
an earlier one as *the model* does not make the earlier fields disposable,
because the comparison between them is a finding.

### SLIM — archive the recipe, delete the fields

| case | size at audit | superseded by |
|---|---|---|
| `waam-arc-01-0p5mm-790W` | 5.1 G | comparison retracted as unusable (4 disqualifying problems) |
| `waam-arc-01-A2-arcpressure` | 1.8 G | P1114W; text already in `results/A2-arcpressure-790W` |
| `waam-arc-01-A1ref` | 1.1 G | A2, then P1114W |
| `waam-arc-01-half-arconly` | 891 M | full-width domain with droplets |
| `waam-arc-01-ref` | 828 M | superseded twice over |
| `waam-arc-01-A1-gz-oldmesh` | 724 M | old mesh, pre gravity-direction fix |
| `waam-arc-01-A1pre-Tclamp` | 569 M | A1ref |
| `waam-arc-01-emON` | 275 M | early EM switch test |
| `waam-arc-01-gminusY-baseline` | 146 M | folded into A1 |
| `waam-arc-01-3mm-argon` | 53 M | early, 2 writes |
| `Plate2D-2b-dt2` | 92 M | timestep check, conclusion banked |
| `Plate2D-mhd-2branch` | 48 M | two-branch kappa dev, superseded by postpatch |
| `Plate2D-mhd-noise` | 44 M | noise-floor study; the 25% `max|U|` floor is banked |
| `Plate2D` | 31 M | pre-MHD dev case |
| `Plate2D-nomhd` | 27 M | superseded by the Gate B pair |
| `P2D-par` / `P2D-ser2` / `P2D-ser` | 63 M | decomposition checks |
| `hartmann1d-fine` | 2.9 M | mesh-refinement check, conclusion banked |

≈ **11.7 GB**.

### PURGE — whole directory

`zz-lorentz-smoke` (Lorentz patch smoke test), `waam-arc-01-A1ref-1512W`
(116 K, 0 time directories, never ran).

---

## 4. Procedure

### Phase 0 — before touching anything

No solver may be running, and every reference case must have its comparison
text banked in `~/thesis/results/`. Do not clean a reference case until its
successor has finished and been compared.

Also re-extract any derived text still living in `/tmp`. The compaction step
in Phase 3 needs `wsl --shutdown`, and `/tmp` is not guaranteed to survive a
restart.

```bash
pgrep -x laserbeamFoam    # must print nothing
```

### Phase 1 — measure, then dedupe

A parallel case stores every write twice: once in `processor*/<t>`, and again
at the case root after `reconstructPar`. The parallel copy is redundant **if
and only if** the reconstructed copy exists, so check per time directory,
never in bulk. `processor*/0` and `processor*/constant` stay, so the
decomposition remains valid, and the parallel form is one
`decomposePar -fields -time <t>` away if ever needed.

**M3 and M4 exist almost entirely in parallel form.** Only 1.06–1.20 were
reconstructed. The per-time check below protects them, but do not dedupe
either one expecting large savings.

```bash
cd ~/thesis/cases/00-tutorials
for c in waam-arc-01-P1114W waam-arc-01-M1-mushy waam-arc-01-M2-full waam-arc-01-M3-darcyA1e5 waam-arc-01-M4-1512W; do
  [ -d "$c/processor0" ] || continue
  dup=0; orph=""
  for t in $(ls -d "$c"/processor0/[0-9]* 2>/dev/null | xargs -n1 basename | grep -v orig); do
    [ "$t" = "0" ] && continue
    if [ -d "$c/$t" ]; then dup=$((dup+1)); else orph="$orph $t"; fi
  done
  echo "$c: $dup dup, proc $(du -shc $c/processor* 2>/dev/null|tail -1|cut -f1), parallel-only:${orph:- none}"
done
```

`parallel-only: none` means the whole `processor*` tree is a second copy and
can go without losing anything. Anything listed there exists **only** in
parallel form and must be reconstructed before its processor copy is touched.

Then, for each case that reported `none` (set `c` to the case name first):

```bash
for t in $(ls -d "$c"/processor0/[0-9]* | xargs -n1 basename | grep -v orig); do
  [ "$t" = "0" ] && continue
  [ -d "$c/$t" ] && rm -rf "$c"/processor*/"$t"
done
```

### Phase 2 — slim the superseded runs

This archives the recipe to `~/thesis/results/archive/<case>/` first, verifies
the archive, refuses anything protected, and refuses any case with a live
solver in its cwd. Logs are kept as head+tail excerpts rather than gzipped in
full, so the archive step itself cannot consume meaningful space.

```bash
cd ~/thesis/cases/00-tutorials || exit 1
KEEP=" waam-arc-01-P1114W waam-arc-01-M1-mushy waam-arc-01-M2-full waam-arc-01-M3-darcyA1e5 waam-arc-01-M4-1512W hartmann1d-postpatch hartmann1d Plate2D-mhd-postpatch Plate2D-mhd "
SLIM="waam-arc-01-0p5mm-790W waam-arc-01-A2-arcpressure waam-arc-01-A1ref
waam-arc-01-half-arconly waam-arc-01-ref waam-arc-01-A1-gz-oldmesh
waam-arc-01-A1pre-Tclamp waam-arc-01-emON waam-arc-01-gminusY-baseline
waam-arc-01-3mm-argon Plate2D-2b-dt2 Plate2D-mhd-2branch Plate2D-mhd-noise
Plate2D Plate2D-nomhd P2D-par P2D-ser2 P2D-ser hartmann1d-fine
zz-lorentz-smoke waam-arc-01-A1ref-1512W"
LIVE=$(for p in $(pgrep -x laserbeamFoam); do readlink -f /proc/$p/cwd 2>/dev/null; done | sort -u)
A=~/thesis/results/archive; mkdir -p "$A"
for c in $SLIM; do
  [ -d "$c" ] || continue
  case "$KEEP" in *" $c "*) echo "REFUSE protected: $c"; continue;; esac
  case "$LIVE" in *"$PWD/$c"*) echo "REFUSE live: $c"; continue;; esac
  d="$A/$c"; mkdir -p "$d"
  tar -C "$c" --exclude=polyMesh -cf - system constant 0.orig 2>/dev/null | tar -C "$d" -xf - 2>/dev/null
  for f in "$c"/Allrun "$c"/Allclean "$c"/*.sh; do [ -f "$f" ] && cp -a "$f" "$d/"; done
  for L in "$c"/log*; do [ -f "$L" ] && { head -c 200000 "$L"; echo; echo "...ELIDED..."; tail -c 400000 "$L"; } > "$d/$(basename "$L").headtail"; done
  { echo "case     $c"; echo "archived $(date)"; du -sh "$c"; echo "times:";
    ls -d "$c"/[0-9]* 2>/dev/null | grep -v orig | xargs -r -n1 basename; } > "$d/MANIFEST.txt"
  [ -s "$d/MANIFEST.txt" ] || { echo "ARCHIVE FAILED, skipping $c"; continue; }
  b=$(du -sk "$c" | cut -f1)
  ls -d "$c"/[0-9]* 2>/dev/null | grep -v orig | xargs -r rm -rf
  rm -rf "$c"/processor* "$c"/constant/polyMesh "$c"/postProcessing "$c"/log*
  echo "slimmed $c   $(( (b - $(du -sk "$c"|cut -f1)) / 1024 ))M freed"
done
echo; df -h ~ | tail -1
```

### Phase 3 — return the bytes to Windows

Nothing above has returned a single byte to Windows yet. Three steps, in order.

**3a. Discard the freed blocks (inside WSL).** ext4 does not zero freed
blocks, so the vhdx has no idea they are free. `fstrim` issues the discards
that punch the holes. Without this, compaction reclaims almost nothing.

```bash
sudo fstrim -av
```

**3b. Stop the distro (Windows PowerShell).** This cannot be done from inside
Ubuntu. Close every Ubuntu terminal and VS Code WSL window first.

```powershell
wsl --shutdown
wsl --list --verbose     # must show Stopped
```

**3c. Compact with diskpart (PowerShell run as administrator).**

**Do not use `wsl --manage Ubuntu-22.04 --set-sparse true`.** Recent WSL
builds refuse it unless `--allow-unsafe` is added, because sparse VHDs carry a
known data-corruption risk. This vhdx holds the only copy of the case fields,
so never force it.

```
diskpart
  select vdisk file="D:\WSL\Ubuntu-22.04\ext4.vhdx"
  attach vdisk readonly
  compact vdisk
  detach vdisk
  exit
```

**Verify (PowerShell):**

```powershell
'{0:N1} GB' -f ((Get-Item 'D:\WSL\Ubuntu-22.04\ext4.vhdx').Length/1GB)
Get-PSDrive C, D | Format-Table Name, @{n='FreeGB';e={[math]::Round($_.Free/1GB,1)}} -Auto
```

Give each command its own explicit `Format-Table`. When several commands are
pasted together without it, PowerShell reuses the first table's columns for
the later output and prints **blank rows**. That happened during the D: move
check.

---

## 5. Host-side levers, for reference

Not WSL, but this is where the rest of C:'s 200 GB went.

| target | size | recovers | note |
|---|---|---|---|
| `pagefile.sys` | 19.8 G | ~12 G | cap at 8–12 GB, needs reboot. Do **not** cut to minimum: 6 MPI ranks can spike |
| `AppData\Local\Arduino15` | 7.9 G | 5–7 G | stale board cores, unrelated to the thesis |
| `anaconda3` | 11.6 G | 2–5 G | `conda clean --all -y` removes only the package-tarball cache and cannot break an env |
| `AppData\Local\Packages` | 25.4 G | — | 9.0 G of this is the Claude desktop sandbox (`rootfs.vhdx`); leave it while sessions are linked |
| `C:\Windows` (WinSxS) | 25.8 G | 3–8 G | `Dism /Online /Cleanup-Image /StartComponentCleanup /ResetBase`, admin |

Hibernation is already disabled, so `powercfg -h off` is not available as a
lever here.

---

## 6. Recovering a slimmed case

Set `c` to the case name, then:

```bash
cp -r ~/thesis/results/archive/"$c" ~/thesis/cases/00-tutorials/"$c"
cd ~/thesis/cases/00-tutorials/"$c"
cp -r 0.orig 0 && blockMesh && setFields
decomposePar && mpirun -np 6 laserbeamFoam -parallel
```

`MANIFEST.txt` in each archive lists the time directories that existed and the
size the case was, so a re-run can be checked against the original. The recipe
is also on GitHub at tag `0p5mm-M1-M4` for every case that existed then.
