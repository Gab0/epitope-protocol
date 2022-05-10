import re


def process_simulation_name(name: str) -> str:
    """
    Convert internal simulation codes into readable names.
    Useful for chart and table labels or whatever else.
    """

    # Parse single mutation name patterns.
    mutation_pat = re.findall(r"mutation_(\w\d+\w)", name)

    if mutation_pat:
        mutation = mutation_pat[0]
        return f"Mutação {mutation}"

    # Parse variation name patterns.
    number_pat = re.findall(r"\d+-{0,1}\d*", name)
    identifier_pat = re.findall(r"^[^\d]+", name)

    if number_pat:
        number = number_pat[0]
    else:
        number = ""

    domain_pattern = r"^(\w)_"
    domains = re.findall(domain_pattern, name)

    if domains:
        name += f", domínio {domains[0]}"
        name = re.sub(domain_pattern, "", name)

    if identifier_pat:
        identifier_map = {
            "DUMMY": "Artificial",
            "RANDOM": "Aleatória",
            "NAT": "Natural",
            "mutate": "Natural",
            "mutation": "Natural",
            "MUTATE": "Natural",
            "DERIV": "Derivada"
        }
        identifier_code = identifier_pat[0].split("_")[-1]
        try:
            identifier = identifier_map[identifier_code]
        except KeyError:
            identifier = "Desconhecido"

    else:
        return name

    if number == "0":
        return "Original"

    return " ".join([identifier, number])
