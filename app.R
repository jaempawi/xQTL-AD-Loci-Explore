library(shiny)
library(DT)
library(dplyr)
library(readr)
library(stringr)
library(bslib)

# ─── Load data ───────────────────────────────────────────────────────────────
dat <- read_csv("data.csv", show_col_types = FALSE) %>%
  mutate(
    chr = as.integer(chr),
    significance = factor(significance, levels = c("genome wide", "suggestive", "ns")),
    top_confidence = factor(top_confidence, levels = c("CL1","CL2","CL3","CL4","CL5"))
  )

# Trans columns — present after pipeline is re-run with trans data
has_trans <- all(c("n_trans_genes","trans_genes") %in% names(dat))
if (!has_trans) {
  dat$n_trans_genes   <- NA_integer_
  dat$trans_genes     <- NA_character_
  dat$n_trans_contexts <- NA_integer_
  dat$trans_contexts  <- NA_character_
}

ct_cols <- c("ct_Brain_xQTL","ct_Exc_xQTL","ct_Inh_xQTL","ct_Oli_xQTL",
             "ct_OPC_xQTL","ct_Ast_xQTL","ct_Microglia_xQTL","ct_Bulk_Immune_xQTL")
ct_labels <- c("Brain","Excitatory","Inhibitory","Oligodendrocyte",
               "OPC","Astrocyte","Microglia","Bulk Immune")
ct_colors <- c("#7c3aed","#0891b2","#9333ea","#059669",
               "#0d9488","#d97706","#dc2626","#be185d")

gwas_studies <- sort(unique(unlist(strsplit(na.omit(dat$gwas_assoc), "\\|"))))
gwas_studies <- gwas_studies[nchar(trimws(gwas_studies)) > 0]

# ─── Confidence color helper ─────────────────────────────────────────────────
conf_badge <- function(x) {
  cols <- c(CL1="#052e16", CL2="#14532d", CL3="#15803d", CL4="#16a34a", CL5="#d1fae5")
  tcols <- c(CL1="#4ade80", CL2="#86efac", CL3="#bbf7d0", CL4="#dcfce7", CL5="#065f46")
  ifelse(is.na(x), "—",
    sprintf('<span style="background:%s;color:%s;padding:2px 7px;border-radius:4px;font-size:10px;font-weight:700;font-family:monospace">%s</span>',
      cols[x], tcols[x], x))
}

sig_badge <- function(x) {
  bg <- c("genome wide"="#fefce8","suggestive"="#fff7ed","ns"="#f8f9fa")
  tc <- c("genome wide"="#713f12","suggestive"="#9a3412","ns"="#6c757d")
  bd <- c("genome wide"="#fde047","suggestive"="#fb923c","ns"="#dee2e6")
  lbl <- c("genome wide"="GW","suggestive"="Sug","ns"="NS")
  ifelse(is.na(x), "—",
    sprintf('<span style="background:%s;color:%s;border:1px solid %s;padding:2px 7px;border-radius:4px;font-size:10px;font-weight:600;font-family:monospace">%s</span>',
      bg[x], tc[x], bd[x], lbl[x]))
}

bool_badge <- function(x) {
  ifelse(is.na(x), '<span style="color:#adb5bd">—</span>',
    ifelse(x == TRUE | x == "TRUE",
      '<span style="color:#198754;font-weight:600">✓</span>',
      '<span style="color:#adb5bd">✗</span>'))
}

dot_html <- function(row) {
  colors <- ct_colors
  names(colors) <- ct_cols
  dots <- sapply(ct_cols, function(c) {
    on <- !is.na(row[[c]]) && row[[c]] == TRUE
    col <- colors[c]
    sprintf('<span style="display:inline-block;width:8px;height:8px;border-radius:50%%;background:%s;opacity:%s;margin-right:2px"></span>',
      col, ifelse(on, "1", "0.12"))
  })
  paste(dots, collapse="")
}

# ─── UI ──────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bg = "#ffffff",
    fg = "#1a1a2e",
    primary = "#2c4a8c",
    base_font = font_google("Lora"),
    code_font = font_google("Source Code Pro"),
    heading_font = font_google("Playfair Display")
  ),

  tags$head(tags$style(HTML("
    /* ── Global ── */
    body { background: #f9f7f4; font-family: 'Lora', Georgia, serif; font-size: 13px; }
    .container-fluid { padding: 0 !important; }

    /* ── Header ── */
    .app-header {
      background: #ffffff;
      border-bottom: 2px solid #e8e0d5;
      padding: 16px 28px;
      display: flex;
      align-items: flex-end;
      justify-content: space-between;
      box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    }
    .app-title { font-family: 'Playfair Display', Georgia, serif; font-size: 20px; font-weight: 700; color: #1a1a2e; letter-spacing: -0.01em; }
    .app-subtitle { font-size: 11px; color: #8a7f72; margin-top: 2px; font-style: italic; font-family: 'Lora', serif; }
    .header-pills { display: flex; gap: 6px; }
    .hpill {
      background: #f0ecf8; border: 1px solid #d6cff0; color: #3d2b8c;
      border-radius: 20px; padding: 4px 12px; font-size: 11px;
      font-family: 'Source Code Pro', monospace; font-weight: 500;
    }
    .hpill strong { color: #2c4a8c; }

    /* ── Layout ── */
    .layout-wrap { display: flex; height: calc(100vh - 70px); }

    /* ── Sidebar ── */
    .sidebar-panel {
      width: 260px; min-width: 260px;
      background: #ffffff;
      border-right: 1px solid #e8e0d5;
      overflow-y: auto;
      padding: 18px 16px;
      display: flex; flex-direction: column; gap: 18px;
    }
    .filter-section { display: flex; flex-direction: column; gap: 8px; }
    .filter-title {
      font-family: 'Source Code Pro', monospace;
      font-size: 9px; font-weight: 600;
      text-transform: uppercase; letter-spacing: .1em;
      color: #8a7f72;
      padding-bottom: 5px; border-bottom: 1px solid #e8e0d5;
    }
    .form-control, .form-select {
      font-family: 'Source Code Pro', monospace !important;
      font-size: 11px !important;
      border: 1.5px solid #e8e0d5 !important;
      background: #faf8f5 !important;
      color: #1a1a2e !important;
      border-radius: 6px !important;
    }
    .form-control:focus, .form-select:focus {
      border-color: #2c4a8c !important;
      box-shadow: 0 0 0 3px rgba(44,74,140,0.1) !important;
    }

    /* Checkboxes */
    .shiny-input-checkboxgroup label { font-size: 11px !important; color: #4a4060 !important; }
    .shiny-input-checkboxgroup .checkbox { margin: 2px 0 !important; }

    /* Reset button */
    #reset_btn {
      background: transparent; border: 1.5px solid #e8e0d5; color: #8a7f72;
      font-family: 'Source Code Pro', monospace; font-size: 11px;
      border-radius: 6px; width: 100%; padding: 6px;
      cursor: pointer; transition: all .15s;
    }
    #reset_btn:hover { border-color: #c53030; color: #c53030; background: #fff5f5; }

    /* ── Main panel ── */
    .main-panel { flex: 1; overflow: hidden; display: flex; flex-direction: column; }

    /* ── Toolbar ── */
    .toolbar {
      background: #ffffff; border-bottom: 1px solid #e8e0d5;
      padding: 10px 18px; display: flex; align-items: center;
      justify-content: space-between; gap: 10px; flex-shrink: 0;
    }
    .toolbar-left { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
    .result-count { font-family: 'Source Code Pro', monospace; font-size: 11px; color: #8a7f72; }
    .result-count strong { color: #2c4a8c; }
    .col-toggle-label { font-family: 'Source Code Pro', monospace; font-size: 10px; color: #adb5bd; text-transform: uppercase; letter-spacing: .06em; }
    .ct-toggle-btn {
      padding: 3px 9px; border-radius: 4px; font-size: 10px;
      font-family: 'Source Code Pro', monospace; cursor: pointer;
      border: 1.5px solid #e8e0d5; background: transparent;
      color: #8a7f72; transition: all .15s; font-weight: 500; margin: 1px;
    }
    .ct-toggle-btn.active { color: #fff; border-color: transparent; }
    .ct-toggle-btn.active.brain { background: #7c3aed; }
    .ct-toggle-btn.active.exc   { background: #0891b2; }
    .ct-toggle-btn.active.inh   { background: #9333ea; }
    .ct-toggle-btn.active.oli   { background: #059669; }
    .ct-toggle-btn.active.opc   { background: #0d9488; }
    .ct-toggle-btn.active.ast   { background: #d97706; }
    .ct-toggle-btn.active.mic   { background: #dc2626; }
    .ct-toggle-btn.active.imm   { background: #be185d; }
    #download_csv {
      background: #2c4a8c; border: none; color: #fff; padding: 6px 14px;
      font-family: 'Source Code Pro', monospace; font-size: 11px; font-weight: 600;
      border-radius: 6px; cursor: pointer;
      box-shadow: 0 2px 6px rgba(44,74,140,0.25);
    }
    #download_csv:hover { background: #1e3a6e; }

    /* ── Table ── */
    .table-wrap { flex: 1; overflow: auto; padding: 0; }
    table.dataTable { font-family: 'Source Code Pro', monospace !important; font-size: 11px !important; border-collapse: collapse !important; }
    table.dataTable thead th {
      background: #f4f1eb !important; color: #4a4060 !important;
      font-family: 'Source Code Pro', monospace !important; font-size: 10px !important;
      font-weight: 600 !important; text-transform: uppercase !important;
      letter-spacing: .06em !important; border-bottom: 2px solid #d6cff0 !important;
      border-right: 1px solid #e8e0d5 !important; padding: 8px 10px !important;
      white-space: nowrap !important;
    }
    table.dataTable tbody td {
      border-right: 1px solid #f0ecf8 !important; padding: 6px 10px !important;
      border-bottom: 1px solid #f4f1eb !important; vertical-align: middle !important;
    }
    table.dataTable tbody tr:hover td { background: #f8f5ff !important; }
    table.dataTable tbody tr.selected td { background: #ede9fe !important; }
    .dataTables_wrapper .dataTables_info,
    .dataTables_wrapper .dataTables_length,
    .dataTables_wrapper .dataTables_filter,
    .dataTables_wrapper .dataTables_paginate {
      font-family: 'Source Code Pro', monospace !important; font-size: 11px !important; color: #8a7f72 !important;
    }
    .dataTables_wrapper .dataTables_paginate .paginate_button.current {
      background: #2c4a8c !important; color: #fff !important; border-color: #2c4a8c !important;
      border-radius: 4px !important;
    }
    .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
      background: #f0ecf8 !important; color: #2c4a8c !important; border-color: #d6cff0 !important;
    }
    .dataTables_wrapper { padding: 10px 16px !important; }
    .dataTables_wrapper .dataTables_filter input {
      font-family: 'Source Code Pro', monospace !important; font-size: 11px !important;
      border: 1.5px solid #e8e0d5 !important; border-radius: 6px !important;
      padding: 4px 8px !important;
    }

    /* ── Detail panel ── */
    .detail-panel {
      position: fixed; right: 0; top: 0; bottom: 0; width: 430px;
      background: #ffffff; border-left: 2px solid #e8e0d5;
      z-index: 1000; overflow-y: auto;
      transform: translateX(100%); transition: transform .28s cubic-bezier(.4,0,.2,1);
      box-shadow: -8px 0 32px rgba(0,0,0,0.08);
    }
    .detail-panel.open { transform: translateX(0); }
    .detail-header {
      padding: 18px 20px; border-bottom: 1px solid #e8e0d5;
      background: linear-gradient(to bottom, #f4f1eb, #ffffff);
      position: sticky; top: 0; z-index: 1;
    }
    .detail-gene { font-family: 'Playfair Display', serif; font-size: 22px; font-weight: 700; color: #2c4a8c; }
    .detail-variant { font-family: 'Source Code Pro', monospace; font-size: 10px; color: #8a7f72; margin-top: 4px; line-height: 1.6; }
    .detail-close {
      position: absolute; top: 14px; right: 16px;
      background: transparent; border: 1.5px solid #e8e0d5; color: #8a7f72;
      width: 28px; height: 28px; border-radius: 50%; cursor: pointer; font-size: 14px;
      display: flex; align-items: center; justify-content: center; transition: all .15s;
    }
    .detail-close:hover { border-color: #c53030; color: #c53030; }
    .detail-body { padding: 16px 20px; }
    .detail-section { margin-bottom: 18px; }
    .detail-section-title {
      font-family: 'Source Code Pro', monospace; font-size: 9px; font-weight: 600;
      text-transform: uppercase; letter-spacing: .1em; color: #8a7f72;
      padding-bottom: 6px; border-bottom: 1px solid #e8e0d5; margin-bottom: 10px;
    }
    .detail-row { display: flex; justify-content: space-between; padding: 4px 0; gap: 8px; border-bottom: 1px solid #faf8f5; }
    .detail-key { font-size: 11px; color: #8a7f72; min-width: 140px; font-family: 'Lora', serif; }
    .detail-val { font-family: 'Source Code Pro', monospace; font-size: 11px; color: #1a1a2e; text-align: right; word-break: break-all; }
    .ctx-item {
      padding: 6px 10px; border-radius: 6px; font-family: 'Source Code Pro', monospace;
      font-size: 10px; border: 1.5px solid #e8e0d5; color: #4a4060;
      background: #faf8f5; margin-bottom: 4px; line-height: 1.5;
    }
    .ctx-item.cl1 { border-color: #166534; background: #f0fdf4; color: #052e16; }
    .ctx-item.cl2 { border-color: #16a34a; background: #f0fdf4; color: #14532d; }
    .ctx-item.cl3 { border-color: #22c55e; background: #f0fdf4; color: #15803d; }
    .ctx-item.cl4 { border-color: #4ade80; background: #f0fdf4; color: #166534; }
    .ctx-item.cl5 { border-color: #86efac; background: #f0fdf4; color: #15803d; }
    .ct-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 5px; margin-top: 4px; }
    .ct-item {
      padding: 5px 9px; border-radius: 5px; font-size: 11px;
      font-family: 'Lora', serif; border: 1.5px solid #e8e0d5;
      display: flex; align-items: center; gap: 5px;
    }
    .ct-item.has { background: #faf8f5; font-weight: 600; color: #1a1a2e; }
    .ct-item.no { color: #c0bdb8; }
    .overlay { display: none; position: fixed; inset: 0; background: rgba(0,0,0,0.15); z-index: 999; backdrop-filter: blur(1px); }
    .overlay.show { display: block; }

    /* Scrollbar */
    ::-webkit-scrollbar { width: 5px; height: 5px; }
    ::-webkit-scrollbar-track { background: transparent; }
    ::-webkit-scrollbar-thumb { background: #d6cff0; border-radius: 3px; }
    ::-webkit-scrollbar-thumb:hover { background: #2c4a8c; }
  "))),

  # Header
  div(class="app-header",
    div(
      div(class="app-title", "xQTL AD Loci Explorer"),
      div(class="app-subtitle", "FunGen Consortium · ADSP Functional Genomics · Alzheimer's Disease Sequencing Project")
    ),
    div(class="header-pills",
      div(class="hpill", "Loci ", tags$strong(textOutput("stat_loci", inline=TRUE))),
      div(class="hpill", "Variants ", tags$strong(textOutput("stat_vars", inline=TRUE))),
      div(class="hpill", "Genes ", tags$strong(textOutput("stat_genes", inline=TRUE))),
      div(class="hpill", "Filtered ", tags$strong(textOutput("stat_shown", inline=TRUE)))
    )
  ),

  div(class="layout-wrap",

    # Sidebar
    div(class="sidebar-panel",

      div(class="filter-section",
        div(class="filter-title", "Search"),
        textInput("q", NULL, placeholder="Gene · rsID · locus · GWAS…", width="100%")
      ),

      div(class="filter-section",
        div(class="filter-title", "Chromosome"),
        selectInput("fchr", NULL,
          choices=c("All"="", setNames(sort(unique(na.omit(dat$chr))), paste0("chr", sort(unique(na.omit(dat$chr)))))),
          width="100%")
      ),

      div(class="filter-section",
        div(class="filter-title", "GWAS Study"),
        selectInput("fgwas", NULL, choices=c("All studies"="", setNames(gwas_studies, gwas_studies)), width="100%")
      ),

      div(class="filter-section",
        div(class="filter-title", "GWAS Significance"),
        checkboxGroupInput("fsig", NULL,
          choiceNames=list(
            tags$span(HTML('<span style="background:#fefce8;color:#713f12;border:1px solid #fde047;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:600">GW</span> Genome-wide')),
            tags$span(HTML('<span style="background:#fff7ed;color:#9a3412;border:1px solid #fb923c;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:600">Sug</span> Suggestive')),
            tags$span(HTML('<span style="background:#f8f9fa;color:#6c757d;border:1px solid #dee2e6;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:600">NS</span> Not significant'))
          ),
          choiceValues=list("genome wide","suggestive","ns"))
      ),

      div(class="filter-section",
        div(class="filter-title", "Top Confidence Level"),
        checkboxGroupInput("fconf", NULL,
          choiceNames=list(
            tags$span(HTML('<span style="background:#052e16;color:#4ade80;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:700">CL1</span> MR/cTWAS + CS overlap')),
            tags$span(HTML('<span style="background:#14532d;color:#86efac;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:700">CL2</span> MR/cTWAS + coloc')),
            tags$span(HTML('<span style="background:#15803d;color:#bbf7d0;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:700">CL3</span> TWAS + CS/coloc')),
            tags$span(HTML('<span style="background:#16a34a;color:#dcfce7;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:700">CL4</span> CS overlap')),
            tags$span(HTML('<span style="background:#d1fae5;color:#065f46;padding:1px 7px;border-radius:4px;font-family:monospace;font-size:10px;font-weight:700">CL5</span> Any xQTL evidence'))
          ),
          choiceValues=list("CL1","CL2","CL3","CL4","CL5"))
      ),

      div(class="filter-section",
        div(class="filter-title", "Cell Type with xQTL"),
        checkboxGroupInput("fct", NULL,
          choiceNames=as.list(ct_labels),
          choiceValues=as.list(ct_cols))
      ),

      div(class="filter-section",
        div(class="filter-title", "Functional Evidence"),
        checkboxGroupInput("ffe", NULL,
          choiceNames=list("TWAS significant","MR significant","cTWAS significant","Has gene","Has context"),
          choiceValues=list("twas","mr","ctwas","has_gene","has_ctx"))
      ),

      div(class="filter-section",
        div(class="filter-title", "Min Inclusion Score"),
        textInput("fmins", NULL, placeholder="e.g. 0.1", width="100%")
      ),

      actionButton("reset_btn", "↺  Reset all filters", width="100%")
    ),

    # Main panel
    div(class="main-panel",

      # Toolbar
      div(class="toolbar",
        div(class="toolbar-left",
          div(class="result-count", "Showing ", tags$strong(textOutput("rcount", inline=TRUE)), " rows"),
          tags$span(style="width:1px;height:14px;background:#e8e0d5;display:inline-block"),
          div(class="col-toggle-label", "Columns"),
          div(
            lapply(seq_along(ct_cols), function(i) {
              tags$button(
                ct_labels[i],
                class=paste("ct-toggle-btn active", c("brain","exc","inh","oli","opc","ast","mic","imm")[i]),
                id=paste0("toggle_", ct_cols[i]),
                onclick=sprintf("Shiny.setInputValue('ct_toggle', '%s_' + Date.now())", ct_cols[i])
              )
            })
          )
        ),
        downloadButton("download_csv", "⬇  Download CSV")
      ),

      # Table
      div(class="table-wrap",
        DTOutput("main_table", height="100%")
      )
    )
  ),

  # Detail panel
  div(class="overlay", id="overlay", onclick="closePanel()"),
  div(class="detail-panel", id="detail_panel",
    div(class="detail-header",
      tags$button("✕", class="detail-close", onclick="closePanel()"),
      div(class="detail-gene", textOutput("dp_gene")),
      div(class="detail-variant", textOutput("dp_variant"))
    ),
    div(class="detail-body",
      uiOutput("dp_body")
    )
  ),

  # JS
  tags$script(HTML("
    var ctVis = {};
    ['ct_Brain_xQTL','ct_Exc_xQTL','ct_Inh_xQTL','ct_Oli_xQTL','ct_OPC_xQTL','ct_Ast_xQTL','ct_Microglia_xQTL','ct_Bulk_Immune_xQTL'].forEach(function(ct) { ctVis[ct] = true; });

    function closePanel() {
      document.getElementById('detail_panel').classList.remove('open');
      document.getElementById('overlay').classList.remove('show');
    }

    function openPanel() {
      document.getElementById('detail_panel').classList.add('open');
      document.getElementById('overlay').classList.add('show');
    }

    Shiny.addCustomMessageHandler('openDetailPanel', function(msg) { openPanel(); });

    $(document).on('click', '#main_table tbody tr', function() {
      var rowData = [];
      $(this).find('td').each(function() { rowData.push($(this).text()); });
      Shiny.setInputValue('selected_row', $(this).index(), {priority: 'event'});
      openPanel();
    });
  "))
)

# ─── Server ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ct_visible <- reactiveVal(setNames(rep(TRUE, length(ct_cols)), ct_cols))

  observeEvent(input$ct_toggle, {
    ct <- sub("_[0-9]+$", "", input$ct_toggle)
    vis <- ct_visible()
    vis[ct] <- !vis[ct]
    ct_visible(vis)
    # Update button style
    session$sendCustomMessage("toggleCTButton", list(ct=ct, active=vis[ct]))
  })

  observeEvent(input$reset_btn, {
    updateTextInput(session, "q", value="")
    updateSelectInput(session, "fchr", selected="")
    updateSelectInput(session, "fgwas", selected="")
    updateCheckboxGroupInput(session, "fsig", selected=character(0))
    updateCheckboxGroupInput(session, "fconf", selected=character(0))
    updateCheckboxGroupInput(session, "fct", selected=character(0))
    updateCheckboxGroupInput(session, "ffe", selected=character(0))
    updateTextInput(session, "fmins", value="")
  })

  filtered <- reactive({
    d <- dat

    if (nchar(trimws(input$q)) > 0) {
      q <- tolower(trimws(input$q))
      d <- d %>% filter(
        str_detect(tolower(coalesce(gene,"")), q) |
        str_detect(tolower(coalesce(rsid,"")), q) |
        str_detect(tolower(coalesce(ADlocus,"")), q) |
        str_detect(tolower(coalesce(variant_ID,"")), q) |
        str_detect(tolower(coalesce(gwas_assoc,"")), q)
      )
    }
    if (nchar(input$fchr) > 0) d <- d %>% filter(chr == as.integer(input$fchr))
    if (nchar(input$fgwas) > 0) d <- d %>% filter(str_detect(coalesce(gwas_assoc,""), fixed(input$fgwas)))
    if (length(input$fsig) > 0) d <- d %>% filter(significance %in% input$fsig)
    if (length(input$fconf) > 0) d <- d %>% filter(top_confidence %in% input$fconf)
    if (length(input$fct) > 0) {
      for (ct in input$fct) d <- d %>% filter(.data[[ct]] == TRUE)
    }
    if ("twas" %in% input$ffe) d <- d %>% filter(twas_sig == TRUE)
    if ("mr" %in% input$ffe) d <- d %>% filter(mr_sig == TRUE)
    if ("ctwas" %in% input$ffe) d <- d %>% filter(ctwas_sig == TRUE)
    if ("has_gene" %in% input$ffe) d <- d %>% filter(!is.na(gene))
    if ("has_ctx" %in% input$ffe) d <- d %>% filter(!is.na(ordered_contexts) & ordered_contexts != "")

    mins <- suppressWarnings(as.numeric(input$fmins))
    if (!is.na(mins)) d <- d %>% filter(!is.na(max_inclusion) & max_inclusion >= mins)

    d
  })

  # Stats
  output$stat_loci  <- renderText(format(n_distinct(na.omit(filtered()$ADlocus)), big.mark=","))
  output$stat_vars  <- renderText(format(n_distinct(na.omit(filtered()$variant_ID)), big.mark=","))
  output$stat_genes <- renderText(format(n_distinct(na.omit(filtered()$gene)), big.mark=","))
  output$stat_shown <- renderText(format(nrow(filtered()), big.mark=","))
  output$rcount     <- renderText(format(nrow(filtered()), big.mark=","))

  # Table columns
  tbl_data <- reactive({
    d <- filtered()
    vis <- ct_visible()

    # Core display cols
    out <- d %>% transmute(
      Locus       = ADlocus,
      Variant     = variant_ID,
      Gene        = gene,
      rsID        = rsid,
      Chr         = paste0("chr", chr),
      Pos         = format(pos, big.mark=",", scientific=FALSE),
      Significance = sig_badge(as.character(significance)),
      `-log10(p)` = round(log10pval, 2),
      `Max Incl.` = round(max_inclusion, 3),
      `cV2F`      = round(cv2f_score, 3),
      Rank        = paste0("#", cv2f_rank),
      Confidence  = conf_badge(as.character(top_confidence)),
      `#Ctx`      = n_contexts,
      `#Trans`    = ifelse(!is.na(n_trans_genes) & n_trans_genes > 0,
                           as.character(n_trans_genes), "—"),
      `xQTL PIP`  = round(xqtl_max_inclusion, 3),
      `TWAS z`    = round(max_twas_z, 2),
      TWAS        = bool_badge(twas_sig),
      MR          = bool_badge(mr_sig),
      cTWAS       = bool_badge(ctwas_sig),
      `Dist TSS`  = ifelse(!is.na(dist_tss), paste0(ifelse(dist_tss>=0,"+",""), round(dist_tss/1000,1), "kb"), "—"),
      `Cell Types` = sapply(seq_len(nrow(d)), function(i) dot_html(d[i,]))
    )

    # Cell type columns if visible
    for (i in seq_along(ct_cols)) {
      ct <- ct_cols[i]
      if (isTRUE(vis[ct])) {
        col_label <- ct_labels[i]
        col_color <- ct_colors[i]
        out[[col_label]] <- ifelse(
          !is.na(d[[ct]]) & d[[ct]] == TRUE,
          sprintf('<span style="display:inline-block;width:10px;height:10px;border-radius:50%%;background:%s"></span>', col_color),
          '<span style="color:#dee2e6">·</span>'
        )
      }
    }
    out
  })

  output$main_table <- renderDT({
    datatable(
      tbl_data(),
      escape = FALSE,
      selection = "single",
      rownames = FALSE,
      options = list(
        pageLength = 50,
        scrollX = TRUE,
        scrollY = "calc(100vh - 220px)",
        dom = '<"top"lf>rt<"bottom"ip>',
        columnDefs = list(
          list(className="dt-left", targets="_all"),
          list(orderable=FALSE, targets=c(match("Cell Types", names(tbl_data()))-1))
        ),
        language = list(search="", searchPlaceholder="Filter table…"),
        initComplete = JS("function(settings, json) { $(this.api().table().header()).css({'font-family':'Source Code Pro, monospace','font-size':'10px'}); }")
      ),
      class = "display compact hover"
    )
  }, server=TRUE)

  # Download
  output$download_csv <- downloadHandler(
    filename = function() paste0("AD_xQTL_filtered_", Sys.Date(), ".csv"),
    content = function(file) {
      write_csv(filtered() %>% select(ADlocus, variant_ID, gene, rsid, chr, pos,
        significance, log10pval, max_inclusion, cv2f_score, cv2f_rank,
        top_confidence, n_contexts, xqtl_max_inclusion, max_twas_z,
        twas_sig, mr_sig, ctwas_sig, dist_tss, gwas_assoc, ordered_contexts), file)
    }
  )

  # Detail panel
  selected_row_data <- reactive({
    req(input$selected_row)
    idx <- input$selected_row + 1
    d <- filtered()
    if (idx > nrow(d)) return(NULL)
    d[idx,]
  })

  output$dp_gene <- renderText({
    r <- selected_row_data(); if(is.null(r)) return("—")
    if(is.na(r$gene)) "(no gene association)" else r$gene
  })

  output$dp_variant <- renderText({
    r <- selected_row_data(); if(is.null(r)) return("—")
    paste0(coalesce(r$variant_ID,""), "\n", coalesce(r$ADlocus,""), " · ", coalesce(r$rsid,""))
  })

  output$dp_body <- renderUI({
    r <- selected_row_data(); if(is.null(r)) return(NULL)

    dr <- function(k, v) {
      if (is.na(v) || v == "") return(NULL)
      div(class="detail-row",
        span(class="detail-key", k),
        span(class="detail-val", HTML(as.character(v)))
      )
    }

    sig_cls <- switch(as.character(r$significance), "genome wide"="sgw", "suggestive"="ssug", "sns")

    gwas_sec <- div(class="detail-section",
      div(class="detail-section-title", "GWAS Evidence"),
      dr("Locus", r$ADlocus), dr("Variant", r$variant_ID), dr("rsID", r$rsid),
      dr("Chr:Pos", if(!is.na(r$chr)) paste0("chr",r$chr,":",format(r$pos,big.mark=",")) else NA),
      dr("Significance", sig_badge(as.character(r$significance))),
      dr("-log10(p)", if(!is.na(r$log10pval)) round(r$log10pval,3) else NA),
      dr("Max Inclusion Score", if(!is.na(r$max_inclusion)) round(r$max_inclusion,4) else NA),
      dr("cV2F Score", if(!is.na(r$cv2f_score)) round(r$cv2f_score,4) else NA),
      dr("cV2F Rank", if(!is.na(r$cv2f_rank)) paste0("#",r$cv2f_rank) else NA),
      dr("GWAS Studies", r$gwas_assoc)
    )

    gene_sec <- if (!is.na(r$gene)) {
      div(class="detail-section",
        div(class="detail-section-title", "Gene & xQTL Summary"),
        dr("Gene", HTML(paste0('<span style="color:#2c4a8c;font-weight:700;font-size:14px;font-family:Playfair Display,serif">',r$gene,'</span>'))),
        dr("Ensembl ID", r$gene_id),
        dr("Distance from TSS", if(!is.na(r$dist_tss)) paste0(ifelse(r$dist_tss>=0,"+",""),round(r$dist_tss/1000,1),"kb") else NA),
        dr("Top Confidence", conf_badge(as.character(r$top_confidence))),
        dr("# Contexts", r$n_contexts),
        dr("Best Context", r$context), dr("Effect", r$effect),
        dr("xQTL Max Inclusion", if(!is.na(r$xqtl_max_inclusion)) round(r$xqtl_max_inclusion,4) else NA),
        dr("Max TWAS z", if(!is.na(r$max_twas_z)) round(r$max_twas_z,3) else NA),
        dr("TWAS Context", r$max_twas_ctx),
        dr("TWAS Significant", bool_badge(r$twas_sig)),
        dr("MR Significant", bool_badge(r$mr_sig)),
        dr("cTWAS Significant", bool_badge(r$ctwas_sig))
      )
    }

    ctx_sec <- if (!is.na(r$ordered_contexts) && r$ordered_contexts != "") {
      ctxs <- strsplit(r$ordered_contexts, " \\| ")[[1]]
      div(class="detail-section",
        div(class="detail-section-title", "Evidence Contexts — ordered by confidence"),
        tagList(lapply(seq_along(ctxs), function(i) {
          ctx <- ctxs[i]
          m <- regmatches(ctx, regexpr("CL(\\d)", ctx))
          cls <- if(length(m)>0) paste0("cl",substr(m,3,3)) else ""
          div(class=paste("ctx-item", cls),
            HTML(paste0('<span style="opacity:.5;margin-right:6px;font-size:9px">',i,'.</span>', ctx))
          )
        }))
      )
    }

    ct_sec <- div(class="detail-section",
      div(class="detail-section-title", "Cell Type xQTL Evidence"),
      div(class="ct-grid",
        tagList(lapply(seq_along(ct_cols), function(i) {
          ct <- ct_cols[i]; has <- !is.na(r[[ct]]) && r[[ct]] == TRUE
          div(class=paste("ct-item", if(has) "has" else "no"),
            HTML(sprintf('<span style="display:inline-block;width:8px;height:8px;border-radius:50%%;background:%s;opacity:%s;margin-right:5px"></span>',
              ct_colors[i], if(has) "1" else "0.2")),
            ct_labels[i]
          )
        }))
      )
    )

    # Trans xQTL section — shown only when trans data is available
    trans_sec <- if (!is.na(r$n_trans_genes) && r$n_trans_genes > 0) {
      trans_gene_list <- strsplit(r$trans_genes, "\\|")[[1]]
      trans_ctx_list  <- if (!is.na(r$trans_contexts)) strsplit(r$trans_contexts, "\\|")[[1]] else character(0)
      div(class="detail-section",
        div(class="detail-section-title",
          HTML('<span style="color:#7c3aed">&#x21F4;</span> Trans xQTL Target Genes'),
          span(style="font-size:10px;color:#adb5bd;font-weight:400;margin-left:6px",
               paste0(r$n_trans_genes, " gene", if(r$n_trans_genes>1) "s" else "",
                      " · ", ifelse(is.na(r$n_trans_contexts), "?", r$n_trans_contexts),
                      " context", if(!is.na(r$n_trans_contexts) && r$n_trans_contexts>1) "s" else ""))
        ),
        div(style="display:flex;flex-wrap:wrap;gap:4px;margin-top:4px",
          tagList(lapply(trans_gene_list, function(g) {
            span(style="background:#f5f3ff;color:#6d28d9;padding:2px 8px;border-radius:4px;font-size:11px;font-family:monospace;border:1px solid #ddd6fe", g)
          }))
        ),
        if (length(trans_ctx_list) > 0)
          div(style="margin-top:6px;font-size:11px;color:#6b7280",
            strong("Contexts: "), paste(trans_ctx_list, collapse=" · "))
      )
    }

    tagList(gwas_sec, gene_sec, ctx_sec, ct_sec, trans_sec)
  })
}

shinyApp(ui, server)
