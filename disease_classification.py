"""
Classifies each disease in Outbreaks.csv by primary transmission mode, based on
established CDC/WHO disease-page descriptions of how each pathogen is transmitted.

Diseases with multiple plausible routes, an unspecified causative agent, or a mode that
does not map cleanly onto person-to-person geographic spread (e.g. prion disease, isolated
environmental fungal infection) are left as "unclassified" rather than forced into a
category. This is a primary-route classification for stratification purposes, not a claim
that the listed route is the only one by which a disease can spread.
"""

VECTOR_BORNE = {
    "Arthropod-borne viral fever, virus unspecified",
    "Chikungunya mosquito-borne viral fever",
    "Crimean-Congo haemorrhagic fever",
    "Dengue fever, unspecified",
    "Japanese encephalitis virus disease",
    "Leishmaniasis, unspecified",
    "Malaria, unspecified",
    "Mosquito-borne viral encephalitis, unspecified",
    "O'nyong-nyong mosquito-borne viral fever",
    "Oropouche virus disease",
    "Rift Valley fever",
    "St. Louis encephalitis virus infection",
    "Sylvatic yellow fever ",
    "Venezuelan equine fever",
    "West Nile virus infection",
    "Yellow fever, unspecified",
    "Zika virus disease",
    "Epidemic louse-borne typhus fever due to Rickettsia prowazekii",
    "Bubonic plague",
}

RESPIRATORY = {
    " Influenza due to identified zoonotic or pandemic influenza virus",
    "COVID-19",
    "Coronavirus infection, unspecified site",
    "Diphtheria, unspecified",
    "Measles",
    "Meningococcal meningitis",
    "Middle East respiratory syndrome",
    "Severe acute respiratory syndrome",
    "Tuberculosis of the respiratory system",
    "Viral pneumonia",
    "Whooping cough due to Bordetella pertussis",
}

WATERBORNE_FOODBORNE = {
    "Acute hepatitis A",
    "Acute hepatitis E",
    "Acute poliomyelitis, unspecified",
    "Bacterial foodborne intoxications, unspecified",
    "Botulism",
    "Cholera",
    "Dracunculiasis",
    "Escherichia coli",
    "Foodborne staphylococcal intoxication",
    "Infections due to other Salmonella",
    "Intestinal infections due to Shigella",
    "Listeriosis, unspecified",
    "Typhoid fever",
}

# Direct contact, bloodborne/sexual, or zoonotic-environmental routes that do not fit
# the three categories above; kept separate rather than merged into "unclassified"
# because the route IS known, it just is not vector-borne, respiratory, or
# waterborne/foodborne.
OTHER_KNOWN_ROUTE = {
    "Ebola disease": "direct contact with bodily fluids",
    "Marburg disease": "direct contact with bodily fluids",
    "Lassa fever": "contact with rodent excreta; bodily fluids",
    "Human immunodeficiency virus disease without mention of associated disease or condition, clinical stage unspecified": "bloodborne/sexual",
    "Gonococcal infection, unspecified": "sexual",
    "Monkeypox": "direct/close contact",
    "Rabies, unspecified": "animal bite",
    "Anthrax, unspecified": "zoonotic; cutaneous/inhalation/ingestion",
    "Hantavirus pulmonary syndrome": "inhalation of rodent excreta (environmental, not person-to-person)",
    "Haemorrhagic fever with renal syndrome": "rodent-borne (hantavirus family)",
    "Leptospirosis, unspecified": "contact with water/soil contaminated by animal urine",
    "Arenavirus disease, unspecified": "rodent-borne",
    "Tularaemia, unspecified": "multiple routes: tick/deer fly bite, animal contact, inhalation, ingestion",
    "Legionnaires disease": "inhalation of contaminated water aerosols (environmental, not person-to-person)",
}

# Left unclassified: causative agent or route is unspecified in the record, the
# condition is not a transmissible infectious disease in the conventional sense
# (e.g. prion disease, sepsis as a clinical syndrome), or the route is genuinely
# ambiguous across the diseases the ICD-style label could refer to.
UNCLASSIFIED = {
    "Coccidioidomycosis, unspecified",
    "Creutzfeldt-Jakob disease, unspecified",
    "Enteroviral vesicular stomatitis",
    "Enterovirus infection of unspecified site",
    "Infection, unspecified",
    "Infectious gastroenteritis or colitis without specification of infectious agent",
    "Other specified viral diseases",
    "Plague, unspecified  ",
    "Pseudomonas aeruginosa",
    "Sepsis without septic shock",
    "Unspecified viral disease",
    "Unspecified viral haemorrhagic fever",
    "Viral meningitis not elsewhere classified, unspecified",
}


def classify(disease_name):
    if disease_name in VECTOR_BORNE:
        return "vector-borne"
    if disease_name in RESPIRATORY:
        return "respiratory"
    if disease_name in WATERBORNE_FOODBORNE:
        return "waterborne/foodborne"
    if disease_name in OTHER_KNOWN_ROUTE:
        return "other-known-route"
    return "unclassified"


if __name__ == "__main__":
    import pandas as pd

    df = pd.read_csv("Outbreaks.csv")
    diseases = sorted(df["Disease"].dropna().unique())
    counts = {}
    for d in diseases:
        cat = classify(d)
        counts[cat] = counts.get(cat, 0) + 1
        print(f"{cat:24s} {d}")
    print("\nCategory totals (disease count, not record count):")
    for cat, n in sorted(counts.items()):
        print(f"  {cat}: {n}")
