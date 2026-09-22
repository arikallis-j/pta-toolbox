"""Plot sky reconstructions produced by `observations/analyse.py`.

Reads the `.npz` files `analyse.py` writes to `results/` and renders them.
Three CLI modes, mutually exclusive:

- `--freq k`: one figure with four panels for frequency index `k` —
  clean-map Re(psi) on the sky, a dirty-map placeholder (see caveat below),
  the singular-value spectrum, and the power spectrum C_l^rec vs l. Saved to
  `results/map_freq{k}.pdf`.
- `--all`: reads `results/sky_map_all.npz` and produces a grid of clean-map
  panels (one per frequency, capped at 9) saved to
  `results/maps_all_freqs.pdf`, and a C_l^rec(f) heatmap saved to
  `results/power_spectrum_lf.pdf`.
- `--eigenmodes --freq k`: plots the top singular eigenmodes on the sky.

**Caveat — dirty map and eigenmodes need `V_matrix`.** `analyse.py` does
the SVD of the (l,m) design matrix internally but only saves the *singular
values*, not the right singular vectors `V`. Both the dirty map (which lives
in the truncated SVD eigenmode basis and needs `V` to project back to (l,m))
and the eigenmode sky patterns (`scripts/main.py`'s `plot_eigenmodes`, which
plots columns of `V` directly) are unavailable until `analyse.py` is
extended to also save `V` (e.g. behind a `--save-eigenmodes` flag — not
implemented here). Mode 1's dirty-map panel is shown as a placeholder with a
console warning; `--eigenmodes` prints an error and exits if `V_matrix` is
missing from the `.npz`.

**Sky synthesis.** `scripts/main.py` renders sky maps via
`ptatoolbox.utils.transform.inverse_sfft2` (an s2fft inverse transform over
the *full* l = 0..2*nside-1 range in s2fft's own 1D/2D `flm` layout), then
`healpy.mollview` on the resulting pixel array. `analyse.py`'s `psi_lm`
covers only l = 2..l_max in its own flat (l, m) ordering (via
`build_lm_index`, reused here) evaluated with complex `scipy.special.sph_harm_y`
harmonics — not s2fft's layout, and not HEALPix's own real-alm layout that
`healpy.alm2map` expects. So sky maps here are synthesised manually: sum
`Y_lm(pixel) * coeff` over `build_lm_index(l_max)` at pixel centres from
`healpy.pix2ang`, then `healpy.mollview` for the sky projection — same final
visualisation call as `scripts/main.py`, different (but consistent) alm ->
pixel transform, made necessary by `analyse.py`'s basis choice.
"""
from pathlib import Path
from typing import Optional

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import typer
from scipy.special import sph_harm_y as sph_harm

from observations.analyse import build_lm_index, calc_Nl_vec

app = typer.Typer(add_completion=False)


def infer_l_max(n_basis: int) -> int:
    """Invert build_lm_index's mode count N_basis = (l_max+1)**2 - 4."""
    l_max = round((n_basis + 4) ** 0.5 - 1)
    if len(build_lm_index(l_max)) != n_basis:
        raise ValueError(f"N_basis={n_basis} is not consistent with l = 2..l_max for any integer l_max")
    return l_max


def synthesize_map(coeffs, lm_index, nside):
    """Sum coeffs[nu] * Y_lm(pixel) over lm_index at every HEALPix pixel centre."""
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    l_arr = np.array([l for l, _ in lm_index])[None, :]
    m_arr = np.array([m for _, m in lm_index])[None, :]
    Y = sph_harm(l_arr, m_arr, theta[:, None], phi[:, None])
    return Y @ coeffs


def power_spectrum(coeffs, lm_index):
    """C_l^rec = (1/(2l+1)) * sum_m |coeffs|**2, i.e. the mean over each l's m-block."""
    l_arr = np.array([l for l, _ in lm_index])
    ls = np.unique(l_arr)
    cl = np.array([np.mean(np.abs(coeffs[l_arr == l]) ** 2) for l in ls])
    return ls, cl


def _overlay_positions(positions_path: Optional[Path]):
    if positions_path is None:
        return
    positions = np.load(positions_path)
    theta, phi = hp.vec2ang(positions)
    hp.projscatter(theta, phi, marker=".", color="red", s=10)


def plot_single_frequency(result, nside, positions_path, out_path):
    psi_lm = result["psi_lm"]
    dirty_map = result["dirty_map"]
    singular_values = result["singular_values"]
    freq_nhz = float(result["freq"]) * 1e9
    n_modes = int(result["n_modes"])
    fisher_cond_full = float(result["fisher_cond_full"])

    l_max = infer_l_max(len(psi_lm))
    lm_index = build_lm_index(l_max)

    fig = plt.figure(figsize=(12, 10))

    field = synthesize_map(psi_lm, lm_index, nside)
    hp.mollview(field.real, sub=(2, 2, 1), title=f"Clean map Re(Ψ), f = {freq_nhz:.3f} nHz", fig=fig.number)
    hp.graticule()
    _overlay_positions(positions_path)

    print(
        "WARNING: dirty-map sky panel skipped — the dirty map lives in the truncated "
        "SVD eigenmode basis and needs the right singular vectors V to project back to "
        "(l, m), but analyse.py does not currently save V_matrix. See module docstring."
    )
    ax_placeholder = fig.add_subplot(2, 2, 2)
    ax_placeholder.axis("off")
    ax_placeholder.text(
        0.5,
        0.5,
        "Dirty map sky reconstruction skipped:\nneeds V_matrix from analyse.py's SVD\n"
        "(not currently saved — see module docstring)",
        ha="center",
        va="center",
        wrap=True,
        fontsize=10,
    )

    ax_sv = fig.add_subplot(2, 2, 3)
    ax_sv.bar(np.arange(1, len(singular_values) + 1), singular_values)
    if fisher_cond_full > 0:
        noise_level = singular_values[0] / np.sqrt(fisher_cond_full)
        ax_sv.axhline(noise_level, linestyle="--", color="black", label="estimated noise floor")
        ax_sv.legend()
    ax_sv.set_xlabel("mode index")
    ax_sv.set_ylabel("singular value")
    ax_sv.set_title(f"Singular value spectrum (retained rank = {n_modes})")

    ax_cl = fig.add_subplot(2, 2, 4)
    ls, cl = power_spectrum(psi_lm, lm_index)
    ax_cl.bar(ls, cl)
    ax_cl.set_xlabel("l")
    ax_cl.set_ylabel(r"$C_l^{\mathrm{rec}}$")
    ax_cl.set_title("Reconstructed power spectrum")

    plt.tight_layout()
    plt.savefig(out_path)
    plt.show()


def plot_all_frequencies(data, nside, out_dir):
    freqs = data["freqs"]
    psi_lm_all = data["psi_lm"]
    n_freqs = len(freqs)

    l_max = infer_l_max(psi_lm_all.shape[1])
    lm_index = build_lm_index(l_max)

    n_panels = min(n_freqs, 9)
    if n_freqs > 9:
        print(f"WARNING: {n_freqs} frequencies available, showing only the first 9.")

    fig = plt.figure(figsize=(14, 14))
    for i in range(n_panels):
        field = synthesize_map(psi_lm_all[i], lm_index, nside)
        freq_nhz = freqs[i] * 1e9
        hp.mollview(field.real, sub=(3, 3, i + 1), title=f"{freq_nhz:.2f} nHz", fig=fig.number)
    plt.tight_layout()
    plt.savefig(out_dir / "maps_all_freqs.pdf")
    plt.show()

    ls, _ = power_spectrum(psi_lm_all[0], lm_index)
    cl_matrix = np.array([power_spectrum(psi_lm_all[i], lm_index)[1] for i in range(n_freqs)]).T
    freqs_nhz = freqs * 1e9

    fig2, ax = plt.subplots(figsize=(10, 6))
    mesh = ax.pcolormesh(freqs_nhz, ls, np.log10(cl_matrix), shading="auto", cmap="viridis")
    fig2.colorbar(mesh, ax=ax, label=r"$\log_{10} C_l^{\mathrm{rec}}$")
    ax.set_xlabel("frequency [nHz]")
    ax.set_ylabel("l")
    ax.set_title("Reconstructed power spectrum vs frequency")
    plt.tight_layout()
    plt.savefig(out_dir / "power_spectrum_lf.pdf")
    plt.show()


def plot_eigenmodes(result, nside, n, out_path):
    if "V_matrix" not in result.files:
        print(
            "ERROR: this .npz does not contain V_matrix (the SVD right singular vectors). "
            "analyse.py does not currently save it. Re-run analyse.py with a "
            "--save-eigenmodes flag to save V_matrix, then re-run mapping.py --eigenmodes "
            "(that flag does not exist yet and needs to be added to analyse.py)."
        )
        raise typer.Exit(code=1)

    V = result["V_matrix"]
    singular_values = result["singular_values"]
    n = min(n, V.shape[1])

    l_max = infer_l_max(V.shape[0])
    lm_index = build_lm_index(l_max)
    Nl_vec = calc_Nl_vec(lm_index)

    print(f"Singular values: {singular_values[:n]}")

    ncols = 3
    nrows = -(-n // ncols)
    fig = plt.figure(figsize=(5 * ncols, 5 * nrows))
    for k in range(n):
        alm_mode = V[:, k]
        psi_mode = alm_mode * (Nl_vec / 2)
        field = synthesize_map(psi_mode, lm_index, nside)
        hp.mollview(field.real, sub=(nrows, ncols, k + 1), title=f"mode {k+1}, s = {singular_values[k]:.3e}", fig=fig.number)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.show()


def main(
    freq: Optional[int] = typer.Option(None, help="frequency index for mode 1 (or with --eigenmodes)"),
    all_freqs: bool = typer.Option(False, "--all", help="mode 2: all frequencies from sky_map_all.npz"),
    eigenmodes: bool = typer.Option(False, "--eigenmodes", help="mode 3: plot top eigenmodes for --freq k"),
    nside: int = typer.Option(32, help="HEALPix resolution for the synthesised sky maps"),
    n: int = typer.Option(6, help="number of eigenmodes to plot (--eigenmodes)"),
    positions: Optional[Path] = typer.Option(None, help="path to a .npy of pulsar unit vectors, shape (N_psr, 3)"),
    results_dir: Path = typer.Option(Path("results"), "--results-dir", help="directory with analyse.py's outputs"),
):
    if eigenmodes:
        if freq is None:
            raise typer.BadParameter("--eigenmodes requires --freq k")
        result = np.load(results_dir / f"sky_map_freq{freq}.npz")
        plot_eigenmodes(result, nside, n, results_dir / f"eigenmodes_freq{freq}.pdf")
    elif all_freqs:
        data = np.load(results_dir / "sky_map_all.npz")
        plot_all_frequencies(data, nside, results_dir)
    elif freq is not None:
        result = np.load(results_dir / f"sky_map_freq{freq}.npz")
        plot_single_frequency(result, nside, positions, results_dir / f"map_freq{freq}.pdf")
    else:
        raise typer.BadParameter("specify --freq k, --all, or --eigenmodes --freq k")


if __name__ == "__main__":
    typer.run(main)
