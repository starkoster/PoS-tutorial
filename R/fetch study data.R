#' Extract study data from ClinicalTrials.gov
#'
#' @param condition    ClinicalTrials.gov query search term
#' @param endpoint     Search terms for specific endpoint
#' @param study_start  Earliest study start date
#' @param study_end    Latest study completion date
#' @param add_data     Add baseline characteristics and results
#' @param use_parallel Use multicore mclapply to fetch data and extract values
#'
#' @return tibble or list of 2 tibbles (with results and baseline characteristics)
#'
#' @examples
#' # Fetch study data on migraine trials
#' nct = get_studies()

get_studies <- function(condition    = "Migraine",
                        endpoint     = "migraine days|mmd",
                        study_start  = "2010-01-01",
                        study_end    = "2025-12-01",
                        add_data     = TRUE,
                        use_parallel = TRUE){

  # A raw list of studies matching general search terms
  all_studies = ctg_get_fields(
    condition = condition,
    location = NULL,
    intervention = "Drug",
    title = NULL,
    status = c("COMPLETED"),
    page_size = 1000) %>% silence
  # glimpse(all_studies)

  # Fine selection of useful studies
  out_studies = all_studies %>%
    filter(
      `Study Results` == "YES",
      `Start Date` >= study_start,
      `Completion Date` <= study_end,
      grepl("^DRUG\\:*", Interventions)
    ) %>%
    mutate(.outcome_txt = tolower(paste0(`Primary Outcome Measures`, `Secondary Outcome Measures`, `Other Outcome Measures`)))  %>%
    filter(grepl(endpoint, .outcome_txt))

  # fix names
  names(out_studies) = tolower(names(out_studies))
  out_studies = out_studies %>% rename(id = `other ids`, nct=`nct number`)
  names(out_studies) = names(out_studies) %>% str_replace_all("study", "")
  names(out_studies) = names(out_studies) %>% str_replace_all("outcome measures", "")
  names(out_studies) = names(out_studies) %>% str_replace_all("brief", "")
  names(out_studies) = names(out_studies) %>% str_trim()
  names(out_studies) = names(out_studies) %>% str_replace_all(" ", "_")

  if(add_data){
    # Fetch data in json structure
    if(use_parallel){
      idfail = 1:nrow(out_studies)
      tmp = vector(mode = "list", length = nrow(out_studies))
      i=0
      while(length(idfail)>0){
        tmp[idfail] = parallel::mclapply(out_studies$nct[idfail], get_json, mc.cores = parallel::detectCores()-4) %>% silence
        idfail = which((lapply(tmp, length) %>% unlist)==1)
        i = i+1
      }

      # Extract baseline data
      all_bl = parallel::mclapply(tmp, get_baseline,
                                  mc.cores = parallel::detectCores()-4) %>% bind_rows %>% silence

      # Extract results data
      all_res = parallel::mclapply(tmp, get_results,
                                   mc.cores = parallel::detectCores()-4) %>% bind_rows %>% silence

    } else{
      tmp = lapply(out_studies$nct, get_json) %>% silence

      # Extract baseline data
      all_bl = lapply(tmp, get_baseline) %>% bind_rows %>% silence

      # Extract results data
      all_res = lapply(tmp, get_results) %>% bind_rows %>% silence
    }


    out_studies = list(results = all_res, baseline=all_bl)
  }

  return(out_studies)
}


#' Extract study JSON data from ClinicalTrials.gov
#'
#' @param nct       NCT identification code
#' @param maxit     max iterations on retrying data fetch
#' @param verbose   Print study NCT in console
#'
#' @return list of lists JSON structure
#'
#' @examples
#' # Fetch study data on migraine trials
#' js = get_json("NCT01952574")

get_json <- function(nct, maxit=50, verbose=FALSE){

  base_url = "https://clinicaltrials.gov/api/v2/studies/"

    if(verbose)
      message("Getting results for ", nct, "...")

    url  = paste0(base_url, URLencode(nct))
    resp = GET(url, accept_json())

    if (status_code(resp) != 200) {
      warning(sprintf("Failed to fetch %s: HTTP %s", nct, status_code(resp)))
      return(NULL)
    }
    content_text = content(resp, as = "text", encoding = "UTF-8")
    # parse JSON
    i=0; json = 1
    while(length(json)==1 & i<maxit){
      # Data should be a list of length > 1, otherwise there is an error
      json = fromJSON(content_text, simplifyVector = FALSE)
      i = i+1
    }

    json$nct = nct

    return(json)
}

#' Extract baseline data from ClinicalTrials.gov JSON structure
#'
#' @param j       JSON structure (list of lists)
#'
#' @return tibble with baseline characteristics available on ClinicalTrials.gov
#'
#' @examples
#' # Fetch study data on migraine trials
#' js = get_json("NCT01952574")
#' get_baseline(js)

get_baseline <- function(j){
  measures = j$resultsSection$baselineCharacteristicsModule$measures

  # Get baseline values
  blcs = map_dfr(measures, function(m) {

    overview = tibble(nct = j$nct,
                      param = m$title,
                      unit = m$unitOfMeasure,
                      time_frame = m$timeFrame,
                      param_type = m$paramType,
                      limit_type = m$dispersionType)

    # ---- CLASS rows ----
    class_rows = map_dfr(m$classes,
                         function(cls) {
                           map_dfr(cls$categories,
                                   function(cat) {
                                     map_dfr(cat$measurements,
                                             bind_rows)
                                   })
                         })

    bind_cols(overview, class_rows)
  })

  j$resultsSection$baselineCharacteristicsModule$groups %>% bind_rows %>% rename(groupId = id) %>% # descriptions of each arm
    left_join(map_dfr(j$resultsSection$baselineCharacteristicsModule$denoms, function(x) x$counts) %>% rename(n=value)) %>%  # number of patients in each arm
    left_join(blcs) %>% # baseline values
    relocate(nct) %>%
    silence

}

#' Extract study results from ClinicalTrials.gov JSON structure
#'
#' @param j       JSON structure (list of lists)
#'
#' @return tibble with study results available on ClinicalTrials.gov
#'
#' @examples
#' # Fetch study data on migraine trials
#' js = get_json("NCT01952574")
#' get_results(js)

get_results = function(j) {

  # Loop through outcomes
  measures = j$resultsSection$outcomeMeasuresModule$outcomeMeasures
  out = map_dfr(measures, function(m) {
    endpoint = m$title
    unit = m$unitOfMeasure
    time_frame = m$timeFrame
    param_type = m$paramType
    limit_type = m$dispersionType

    # Lookup group names and counts
    groups = setNames(map_chr(m$groups, "title"), map_chr(m$groups, "id"))
    counts = setNames(map_chr(m$denoms[[1]]$counts, "value"), map_chr(m$denoms[[1]]$counts, "groupId"))

    # ---- CLASS rows ----
    class_rows = map_dfr(m$classes, function(cls) {
      map_dfr(cls$categories, function(cat) {
        map_dfr(cat$measurements, function(meas) {
          tibble(
            nct = j$nct,
            scope = "class",
            endpoint = endpoint,
            unit = unit,
            time_frame = time_frame,
            param_type = param_type,
            limit_type = limit_type,
            group_id = meas$groupId,
            group_name = groups[[meas$groupId]],
            n = as.numeric(counts[[meas$groupId]]),
            value = as.numeric(meas$value),
            lower = as.numeric(meas$lowerLimit),
            upper = as.numeric(meas$upperLimit)
          )
        })
      })
    })

    # ---- ANALYSIS rows ----
    analysis_rows = map_dfr(m$analyses, function(a) {
      g_ids = unlist(a$groupIds)  # FIX: unlist groupIds
      tibble(
        nct = j$nct,
        scope = "analysis",
        endpoint = endpoint,
        unit = unit,
        time_frame = time_frame,
        param_type = a$paramType,
        limit_type = "Confidence Interval",
        group_id = paste(g_ids, collapse = "|"),
        group_name = paste(groups[g_ids], collapse = " vs "),
        n = sum(as.numeric(counts[g_ids]), na.rm = TRUE),
        value = as.numeric(a$paramValue),
        lower = as.numeric(a$ciLowerLimit),
        upper = as.numeric(a$ciUpperLimit),
        p_value = a$pValue,
        analysis_method = a$statisticalMethod
      )
    })

    return(bind_rows(class_rows, analysis_rows))
  })

  return(out)
}
