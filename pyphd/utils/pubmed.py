# Utility functions to query pubmed.
# Copyright (C) 2018  Lisa Perus
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# long with this program.  If not, see <https://www.gnu.org/licenses/>.

# System imports
import os

# Third-party imports
import pandas as pd
from Bio import Entrez

# Pyphd imports
from pyphd.constants import SCIMAGOJR_RANKING_DATA


def pubmed_search(query, sort="Best Match", retmax="20"):
    """Returns list of papers IDs corresponding to a query.

    /!\ Code was directly copied from:
        https://marcobonzanini.com/2015/01/12/searching-pubmed-with-python/
    /!\

    Parameters
    ----------
    query: str
        Query for pubmed.
    sort: str
        Type of sorting. E.g : relevance, Best Match.
    retmax: str
        Number of articles retained.

    Returns
    -------
    results: Entrez results
        Results containing (among other things) list of publications IDs
        corresponding to the query.
    """
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort=sort,
                            retmax=retmax,
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results


def fetch_pubmed_articles_details(query, sort="Best Match", retmax="20"):
    """Returns details of papers corresponding to a query.

    /!\ Code was directly copied from:
        https://marcobonzanini.com/2015/01/12/searching-pubmed-with-python/
    /!\

    ex query: '((((fmri) AND alzheimer)) AND
    ("2018/01/01"[Date - Publication] : "3000"[Date - Publication])) '

    Parameters
    ----------
    query: str
        Query for pubmed.
    sort: str
        Type of sorting. E.g : relevance, Best Match.
    retmax: str
        Number of articles retained.

    Returns
    -------
    papers_results: dict
        Dictionnary containing information about papers
        corresponding to the query.
    """
    results = pubmed_search(query, sort, retmax)
    id_list = results['IdList']
    ids = ','.join(id_list)
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    papers_results = Entrez.read(handle)
    return papers_results


def fetch_pubmed_articles_details_by_ids(
        ids_list, sort="Best Match"):
    """Returns details of papers from pubmed id list.

    /!\ Code was directly copied from:
        https://marcobonzanini.com/2015/01/12/searching-pubmed-with-python/
    /!\


    Parameters
    ----------
    ids_list: list of str
        List of pubmed IDs.
    sort: str
        Type of sorting. E.g : relevance, Best Match.
    retmax: str
        Number of articles retained.

    Returns
    -------
    papers_results: dict
        Dictionnary containing information about papers
        corresponding to the query.
    """
    ids = ','.join(ids_list)
    Entrez.email = 'your.email@example.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    papers_results = Entrez.read(handle)
    return papers_results


def sort_articles_by_abstract(papers_info, patterns, logical="and"):
    """Sort articles by looking at specific patterns/keyword in abstracts.

    Parameters
    ----------
    papers_info: dict
        Dictionnary obtained from fetch_pubmed_articles_details, containing
        information about papers corresponding a query.
    patterns: list of str
        patterns to be found.
    logical: str
        logical link between pattern

    Returns
    -------
    papers_to_retain: dict
        Dictionnary containing information about papers
        with corresponding pattern(s) in abstract
    """
    papers = papers_info["PubmedArticle"]
    papers_to_retain = {"PubmedArticle": []}
    for paper in papers:
        retain_paper = False
        found_patterns = {}
        if "Abstract" not in paper["MedlineCitation"]["Article"].keys():
            article_name = paper[
                                "MedlineCitation"]["Article"]["ArticleTitle"]
            print("Could not find abstract for {0}".format(article_name))
            continue
        abstract = paper["MedlineCitation"]["Article"]["Abstract"][
                    "AbstractText"]

        if logical == "and":
            for elt in abstract:
                for pat in patterns:
                    if pat not in found_patterns:
                        found_patterns[pat] = False

                    # > If pattern has already been found, continue
                    if found_patterns[pat]:
                        continue

                    if pat in elt:
                        found_patterns[pat] = True

            # > Check if all patterns have been found
            all_patterns_found_list = []
            for pat, pat_data in found_patterns.items():
                all_patterns_found_list.append(pat_data)
            all_patterns_found = True
            for elt in all_patterns_found_list:
                if not elt:
                    all_patterns_found = False
            if all_patterns_found:
                retain_paper = True

        elif logical == "or":
            found_patterns = False
            for elt in abstract:
                for pat in patterns:
                    if pat in elt:
                        found_patterns = True
            if found_patterns:
                retain_paper = True
        else:
            raise NotImplementedError(
                "More complex logical association than or/and have not been "
                "implemented yet.")

        # Keep article or not
        if retain_paper:
            papers_to_retain["PubmedArticle"].append(paper)

    return papers_to_retain


def sort_articles_by_journal_sjr(papers_info, outdir):
    """Sort articles by journal sjr.

    Parameters
    ----------
    papers_info: dict
        Dictionnary obtained from fetch_pubmed_articles_details, containing
        information about papers corresponding a query.

    Returns
    -------
    output_csv: str
        Path to file listing articles with their journal rankings.
    outdir: str
        Path to output directory
    """

    # Print info
    print("Using scimagojr file from {0} ...".format(
            SCIMAGOJR_RANKING_DATA["link"]))
    print("Last Download : {0}".format(SCIMAGOJR_RANKING_DATA["date"]))
    print("Version : {0}".format(SCIMAGOJR_RANKING_DATA["version"]))
    papers_ranking = {}

    # Load Scimagorjr data
    scimagojr_ranking_data = pd.read_csv(
        SCIMAGOJR_RANKING_DATA["data"], sep=";")

    # Parse papers
    papers = papers_info["PubmedArticle"]
    for paper in papers:

        # > Get journal ISSN ID
        medline_info = paper["MedlineCitation"]["MedlineJournalInfo"]
        if "ISSNLinking" in medline_info.keys():
            issn_id = medline_info["ISSNLinking"]
        else:
            if "ISSN" in paper["MedlineCitation"]["Article"]["Journal"].keys():
                issn_id = str(
                    paper["MedlineCitation"]["Article"]["Journal"]["ISSN"])
            else:
                issn_id = "Not found"

        # > Scimagorjr does not have an hyphen in its ISSN IDs
        issn_id = issn_id.replace("-", "")

        # > Get journal name
        line_idx = None
        found_journal = 0
        for idx, elt in enumerate(scimagojr_ranking_data["Issn"]):
            if issn_id in elt:
                line_idx = idx

        if line_idx is None:
            if None not in papers_ranking.keys():
                papers_ranking[-1] = {"Name": "No journal found.",
                                      "SJR Best Quartile": "NA",
                                      "H index": "NA",
                                      "Articles": []}
            papers_ranking[-1]["Articles"].append(paper)
        else:
            journal_sjr = scimagojr_ranking_data.iloc[line_idx, :]["SJR"]
            if type(journal_sjr) == str:
                journal_sjr = journal_sjr.replace(",", ".")
            journal_sjr = float(journal_sjr)
            if journal_sjr not in papers_ranking.keys():
                papers_ranking[journal_sjr] = {
                    "Name": scimagojr_ranking_data.iloc[line_idx, :]["Title"],
                    "SJR Best Quartile": scimagojr_ranking_data.iloc[
                        line_idx, :]["SJR Best Quartile"],
                    "H index": scimagojr_ranking_data.iloc[
                        line_idx, :]["H index"],
                    "Articles": []
                }
            papers_ranking[journal_sjr]["Articles"].append(paper)

    # Order ranking and write output
    sorted_ranks = sorted(list(papers_ranking.keys()), reverse=True)

    # Write output file
    output_csv = os.path.join(outdir, "articles_rank_journals.csv")
    if os.path.isfile(output_csv):
        raise ValueError("File '{0}' already exists. Please delete.".format(
            output_csv))
    with open(output_csv, "wt") as open_file:
        open_file.write("Article;Year;Journal;Journal_SJR;")
        open_file.write("Journal_SJR_Best_Quartile;Journal_H_index\n")
        for rank in sorted_ranks:
            for article in papers_ranking[rank]["Articles"]:
                article_name = article[
                    "MedlineCitation"]["Article"]["ArticleTitle"]
                if type(article_name) == list:
                    article_name = article_name[0]
                article_name = article_name.replace(";", " - ")
                article_name = article_name.replace("\"", " ")
                if len(article["MedlineCitation"]["Article"][
                               "ArticleDate"]) == 0:
                    year = "NA"
                else:
                    year = article["MedlineCitation"]["Article"][
                                   "ArticleDate"][0]["Year"]
                line = ";".join([article_name, year,
                                 papers_ranking[rank]["Name"],
                                 str(rank),
                                 papers_ranking[rank]["SJR Best Quartile"],
                                 str(papers_ranking[rank]["H index"])])
                open_file.write(line)
                open_file.write("\n")
    return output_csv
