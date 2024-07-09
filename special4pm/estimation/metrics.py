import math

import mpmath
from numpy import euler_gamma
from scipy.special import digamma

from cachetools import cached


#TODO unify incidence and abundance-based methods in one function
def get_incidence_count(obs_species_counts: dict, i: int) -> int:
    """
    returns the number of species, that have an incidence count of i
    :param obs_species_counts: the species with corresponding incidence counts
    :param i: the incidence count
    :return: the number of species with incidence count i
    """
    return list(obs_species_counts.values()).count(i)


def get_singletons(obs_species_counts: dict) -> int:
    """
    returns the number of singletons species, i.e. those species that have an incidence count of 1
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the number of species with incidence count 1
    """
    return list(obs_species_counts.values()).count(1)


def get_doubletons(obs_species_counts: dict) -> int:
    """
    returns the number of doubleton species, i.e. those species that have an incidence count of 2
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the number of species with incidence count 2
    """
    return list(obs_species_counts.values()).count(2)


def get_number_observed_species(obs_species_counts: dict) -> int:
    """
    returns the number of observed species
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the number of observed species
    """
    return len(obs_species_counts.keys())


def get_total_species_count(obs_species_counts):
    """
    returns the total number of species incidences, i.e. the sum of all species incidences in the reference sample
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the sum of species incidences
    """
    return sum(obs_species_counts.values())


def hill_number(d: int, obs_species_counts: dict) -> float:
    """
    computes sample-based Hill number of order d for the reference sample.
    D=0 - species richness
    D=1 - Exponential of Shannon entropy
    D=2 - Simpson Diversity Index
    :param d: the order of the Hill number
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the sample-based Hill number of order d
    """
    # sample-based species richness
    if d == 0:
        return get_number_observed_species(obs_species_counts)
    # sample-based exponential Shannon diversity
    if d == 1:
        return entropy_exp(obs_species_counts)
    # sample-based Simpson diversity
    if d == 2:
        return simpson_diversity(obs_species_counts)


def entropy_exp(obs_species_counts: dict) -> float:
    """
    computes the exponential of Shannon entropy
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the exponential of Shannon entropy
    """
    total_species_count = get_total_species_count(obs_species_counts)
    # total = sum([obs_species_counts.values])
    return math.exp(-1 * sum(
        [x / total_species_count * math.log(x / total_species_count) for x in obs_species_counts.values()]))


def simpson_diversity(obs_species_counts: dict) -> float:
    """
    computes the Simpson diversity index
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the Simpson diversity index
    """
    total_species_count = get_total_species_count(obs_species_counts)

    a = sum([(x / total_species_count) ** 2 for x in obs_species_counts.values()])
    # TODO check if return 1 is reasonable
    return a ** (1 / (1 - 2)) if a > 0 else 1


'''
Calculate asymptotic Hill number of order d for a reference sample
d=0 Species Richness
d=1 Exponential of Shannon Entropy
d=2 Simpson Diversity Index
'''


def hill_number_asymptotic(d: int, obs_species_counts: dict, sample_size: int, abundance: bool = True) -> float:
    """
    computes asymptotic Hill number of order d for the reference sample, for either abundance data or incidence data.
    D=0 - asymptotic species richness
    D=1 - asymptotic exponential of Shannon entropy
    D=2 - asymptotic Simpson diversity index
    :param d: the order of the Hill number
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :param abundance: flag indicating the data type. Setting this 'True' indicates abundance-based data,
    setting this 'False' indicates incidence-based data
    :return: the asymptotic Hill number of order d
    """
    # asymptotic species richness
    # for species richness, there is no differentiation between abundance and incidence
    if d == 0:
        return estimate_species_richness_chao(obs_species_counts)
    # asymptotic Shannon entropy
    if d == 1:
        if abundance:
            return estimate_exp_shannon_entropy_abundance(obs_species_counts, sample_size)
        # incidence
        else:
            return estimate_exp_shannon_entropy_incidence(obs_species_counts, sample_size)
    # asymptotic Simpson diversity
    if d == 2:
        if abundance:
            return estimate_simpson_diversity_abundance(obs_species_counts, sample_size)
        # incidence
        else:
            return estimate_simpson_diversity_incidence(obs_species_counts, sample_size)


def estimate_species_richness_chao(obs_species_counts: dict) -> float:
    """
    computes the asymptotic(=estimated) species richness using the Chao1 estimator(for abundance data)
    or Chao2 estimator (for incidence data)
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the estimated species richness
    """
    obs_species_count = get_number_observed_species(obs_species_counts)
    f_1 = get_singletons(obs_species_counts)
    f_2 = get_doubletons(obs_species_counts)

    if f_2 != 0:
        return obs_species_count + f_1 ** 2 / (2 * f_2)
    else:
        return obs_species_count + f_1 * (f_1 - 1) / 2


def estimate_species_richness_chao_corrected(obs_species_counts: dict) -> float:
    """
    computes the asymptotic(=estimated) species richness using the Chao1 estimator(for abundance data)
    or Chao2 estimator (for incidence data). Includes a correction term for very small sample sizes
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the estimated species richness
    """
    obs_species_count = get_number_observed_species(obs_species_counts)
    f_1 = get_singletons(obs_species_counts)
    f_2 = get_doubletons(obs_species_counts)

    if f_2 != 0:
        return ((obs_species_count - 1) / obs_species_count) * obs_species_count + f_1 ** 2 / (2 * f_2)
    else:
        return ((obs_species_count - 1) / obs_species_count) * obs_species_count + f_1 * (f_1 - 1) / 2


def estimate_exp_shannon_entropy_abundance(obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the asymptotic(=estimated) exponential of Shannon entropy for abundance-based data
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the estimated exponential of Shannon entropy
    """
    return math.exp(estimate_entropy(obs_species_counts, sample_size))


def estimate_exp_shannon_entropy_incidence(obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the asymptotic(=estimated) exponential of Shannon entropy for incidence-based data
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the estimated exponential of Shannon entropy
    """
    # term h_o is structurally equivalent to abundance based entropy estimation, see eq H7 in appendix H of Hill number paper
    u = sum(obs_species_counts.values())
    h_o = estimate_entropy(obs_species_counts, sample_size)

    return math.exp((sample_size / u) * h_o + math.log(u / sample_size))


def estimate_entropy(obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the estimated Shannon entropy
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the estimated exponential of Shannon entropy
    """
    # TODO make this understandable
    f_1 = get_singletons(obs_species_counts)
    f_2 = get_doubletons(obs_species_counts)

    entropy_known_species = 0

    for x_i in obs_species_counts.values():
        if x_i <= sample_size - 1:
            norm_factor = x_i / sample_size

            #decompose sum(1/x_i,...,1/sample_size) to um(1/1,...,1/sample_size)-sum(1/1,...,1/x_i-1)
            entropy_known_species = entropy_known_species + norm_factor * (harmonic(sample_size) - harmonic(x_i-1))
            #entropy_known_species = entropy_known_species + norm_factor * (mpmath.harmonic(sample_size) - mpmath.harmonic(x_i-1))
            #entropy_known_species = entropy_known_species + norm_factor * sum([1 / k for k in range(x_i, sample_size)])

    a = 0
    if f_2 > 0:
        a = (2 * f_2) / ((sample_size - 1) * f_1 + 2 * f_2)
    if f_2 == 0 and f_1 > 0:
        a = 2 / ((sample_size - 1) * (f_1 - 1) + 2)
    else:
        a = 1

    entropy_unknown_species = 0

    # TODO rethink if this is really necessary
    if a == 1:
        return entropy_known_species

    entropy_unknown_species = (f_1 / sample_size) * ((1 - a) ** (-sample_size + 1)) * (
            -math.log(a) - sum([(1 / r) * ((1 - a) ** r) for r in range(1, sample_size)]))

    return entropy_known_species + entropy_unknown_species


def harmonic(n):
    """Returns an (approximate) value of n-th harmonic number.
    If n>100, use an efficient approximation using the digamma function instead
    http://en.wikipedia.org/wiki/Harmonic_number
    taken from: https://stackoverflow.com/questions/404346/python-program-to-calculate-harmonic-series
     """
    if n <= 100:
        return sum(1/k for k in range(1,n+1))
    else:
        return digamma(n + 1) + euler_gamma


def estimate_simpson_diversity_abundance(obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the asymptotic(=estimated) Simpson diversity for abundance-based data
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the estimated Simpson diversity
    """
    # TODO make this understandable
    denom = 0
    for x_i in obs_species_counts.values():
        if x_i >= 2:
            denom = denom + (x_i * (x_i - 1))
    if denom == 0:
        return 0
    return (sample_size * (sample_size - 1)) / denom


def estimate_simpson_diversity_incidence(obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the asymptotic(=estimated) Simpson diversity for incidence-based data
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the estimated Simpson diversity
    """
    # TODO make this understandable
    u = get_total_species_count(obs_species_counts)
    s = 0

    nom = ((1 - (1 / sample_size)) * u) ** 2

    for y_i in obs_species_counts.values():
        if y_i > 1:
            #    s = s + (sample_size ** 2 * y_i ** 2) / (u ** 2 * sample_size ** 2)
            s = s + (y_i * (y_i - 1))
    if s == 0:
        return 0
    # return s ** (1 / (1 - 2))
    return nom / s


def completeness(obs_species_counts: dict) -> float:
    """
    computes the completeness of the sample data. A value of '1' indicates full completeness,
    whereas as value of '0' indicates total incompleteness
    :param obs_species_counts: the species with corresponding incidence counts
    :return: the estimated completeness
    """
    obs_species_count = get_number_observed_species(obs_species_counts)
    s_P = estimate_species_richness_chao(obs_species_counts)
    if s_P == 0:
        return 0

    return obs_species_count / s_P


def coverage(obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the coverage of the sample data. A value of '1' indicates full coverage,
    whereas as value of '0' indicates no coverage
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the estimated coverage
    """
    f_1 = get_singletons(obs_species_counts)
    f_2 = get_doubletons(obs_species_counts)
    Y = get_total_species_count(obs_species_counts)

    if sample_size == 0:
        return 0
    if f_2 == 0 and sample_size == 1:
        return 0
    if f_1 == 0 and f_2 == 0:
        return 1

    return 1 - f_1 / Y * (((sample_size - 1) * f_1) / ((sample_size - 1) * f_1 + 2 * f_2))


def sampling_effort_abundance(n: float, obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the expected additional sampling effort needed to reach target completeness l for abundance data.
    If f exceeds the current completeness, this function returns 0
    :param n: desired target completeness
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the expected additional sampling effort
    """
    comp = completeness(obs_species_counts)
    f_1 = get_singletons(obs_species_counts)
    f_2 = get_doubletons(obs_species_counts)

    if n <= comp:
        return 0
    if f_2 == 0:
        return 0

    obs_species_count = estimate_species_richness_chao(obs_species_counts)
    #get_number_observed_species(obs_species_counts))

    s_P = 0
    if f_2 != 0:
        s_P = f_1 ** 2 / (2 * f_2)
    else:
        s_P = f_1 * (f_1 - 1) / 2

    return ((sample_size * f_1) / (2 * f_2)) * math.log(s_P / ((1 - n) * (s_P + obs_species_count)))


def sampling_effort_incidence(n: float, obs_species_counts: dict, sample_size: int) -> float:
    """
    computes the expected additional sampling effort needed to reach target completeness l for incidence data.
    If f exceeds the current completeness, this function returns 0
    :param n: desired target completeness
    :param obs_species_counts: the species with corresponding incidence counts
    :param sample_size: the sample size associated with the species incidence counts
    :return: the expected additional sampling effort
    """
    comp = completeness(obs_species_counts)
    f_1 = get_singletons(obs_species_counts)
    f_2 = get_doubletons(obs_species_counts)
    if n <= comp or sample_size < 2:
        return 0
    if f_2 == 0:
        return 0

    # for small sample sizes, correction term is introduced, otherwise math error
    obs_species_count = estimate_species_richness_chao(obs_species_counts)
    #obs_species_count = get_number_observed_species(obs_species_counts)

    s_P = 0
    if f_2 != 0:
        s_P = obs_species_count + (1 - 1 / sample_size) * f_1 ** 2 / (2 * f_2)
    else:
        s_P = obs_species_count + (1 - 1 / sample_size) * f_1 * (f_1 - 1) / 2

    #TODO double check if this is indeed correct

    nom1 = (sample_size / (sample_size - 1))
    nom2 = ((2 * f_2) / (f_1 ** 2))
    nom3 = (n * s_P - obs_species_count)
    nom = (math.log(1 - nom1 * nom2 * nom3))
    denominator = (math.log(1 - ((2 * f_2) / ((sample_size - 1) * f_1 + 2 * f_2))))

    return nom / denominator
