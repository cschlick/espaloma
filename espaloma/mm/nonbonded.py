# =============================================================================
# IMPORTS
# =============================================================================
import torch

# =============================================================================
# CONSTANTS
# =============================================================================
import espaloma as esp

# CODATA 2018
# ref https://en.wikipedia.org/wiki/Coulomb_constant
# K_E = (
#     8.9875517923 * 1e9
#     * unit.kilogram
#     * (unit.meter ** 3)
#     * (unit.second ** (-4))
#     * (unit.angstrom ** (-2))
#     ).value_in_unit(esp.units.COULOMB_CONSTANT_UNIT)


# =============================================================================
# UTILITY FUNCTIONS FOR COMBINATION RULES FOR NONBONDED
# =============================================================================
def geometric_mean(msg="m", out="epsilon"):
    def _geometric_mean(nodes):
        return {out: torch.prod(nodes.mailbox[msg], dim=1).pow(0.5)}

    return _geometric_mean


def arithmetic_mean(msg="m", out="sigma"):
    def _arithmetic_mean(nodes):
        return {out: torch.sum(nodes.mailbox[msg], dim=1).mul(0.5)}

    return _arithmetic_mean


# =============================================================================
# COMBINATION RULES FOR NONBONDED
# =============================================================================
def lorentz_berthelot(g, suffix=""):
    import dgl
    g.multi_update_all(
        {
            "n1_as_%s_in_%s"
            % (pos_idx, term): (
                dgl.function.copy_src(
                    src="epsilon%s" % suffix, out="m_epsilon"
                ),
                geometric_mean(msg="m_epsilon", out="epsilon%s" % suffix),
            )
            for pos_idx in [0, 1]
            for term in ["nonbonded", "onefour"]
        },
        cross_reducer="sum",
    )

    g.multi_update_all(
        {
            "n1_as_%s_in_%s"
            % (pos_idx, term): (
                dgl.function.copy_src(src="sigma%s" % suffix, out="m_sigma"),
                arithmetic_mean(msg="m_sigma", out="sigma%s" % suffix),
            )
            for pos_idx in [0, 1]
            for term in ["nonbonded", "onefour"]
        },
        cross_reducer="sum",
    )

    return g


# =============================================================================
# ENERGY FUNCTIONS
# =============================================================================
def lj_12_6(x, sigma, epsilon):
    """Lennard-Jones 12-6.

    Parameters
    ----------
    x : `torch.Tensor`, `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    sigma : `torch.Tensor`, `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    epsilon : `torch.Tensor`,
        `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    Returns
    -------
    u : `torch.Tensor`, `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    """

    return esp.mm.functional.lj(x=x, sigma=sigma, epsilon=epsilon)


def lj_9_6(x, sigma, epsilon):
    """Lennard-Jones 9-6.

    Parameters
    ----------
    x : `torch.Tensor`, `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    sigma : `torch.Tensor`, `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    epsilon : `torch.Tensor`,
        `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    Returns
    -------
    u : `torch.Tensor`, `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`

    """

    return esp.mm.functional.lj(
        x=x, sigma=sigma, epsilon=epsilon, order=[9, 6], coefficients=[2, 3]
    )


#
# def columb(x, q_prod, k_e=K_E):
#     """ Columb interaction without cutoff.
#
#     Parameters
#     ----------
#     x : `torch.Tensor`, shape=`(batch_size, 1)` or `(batch_size, batch_size, 1)`
#     q_prod : `torch.Tensor`,
#         `shape=(batch_size, 1) or `(batch_size, batch_size, 1)`
#
#     Returns
#     -------
#     u : `torch.Tensor`,
#         `shape=(batch_size, 1)` or `(batch_size, batch_size, 1)`
#
#
#     """
#     return k_e * x / q_prod
