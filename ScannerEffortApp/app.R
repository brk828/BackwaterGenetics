library(shiny)
library(bslib)
library(dplyr)
library(lubridate)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)

# ── Data preparation ──────────────────────────────────────────────────────────

load("data/ReportingData.RData")

month_levels <- c("Jan","Feb","Mar","Apr","May","Jun",
                  "Jul","Aug","Sep","Oct","Nov","Dec")

LocationLabels <- StudyBWAnalysis |>
  count(location, species) |>
  group_by(location) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(LocationLabel = paste0(location, " - ", species)) |>
  select(Location = location, LocationLabel)

ReliableEffort <- StudyBWEffort |>
  mutate(PctDiff = abs(ScanTimeHrs - DeployedHrs) / DeployedHrs) |>
  filter(DeployedHrs > 0, ScanTimeHrs > 0, PctDiff <= 0.15)

SubEffort <- ReliableEffort |>
  rowwise() |>
  mutate(SubMonths = list({
    month_seq <- seq(floor_date(as.Date(Deploy,   tz = "UTC"), "month"),
                     floor_date(as.Date(Retrieve, tz = "UTC"), "month"),
                     by = "month")
    sub_start <- pmax(Deploy,   as.POSIXct(month_seq,                        tz = "UTC"))
    sub_end   <- pmin(Retrieve, as.POSIXct(ceiling_date(month_seq, "month"), tz = "UTC"))
    tibble(
      MonthStart   = as.Date(month_seq),
      SubEffortHrs = as.numeric(difftime(sub_end, sub_start, units = "hours"))
    )
  })) |>
  unnest(SubMonths) |>
  ungroup() |>
  filter(SubEffortHrs > 0) |>
  mutate(
    ScanMonth     = month(MonthStart),
    ScanMonthName = factor(format(MonthStart, "%b"), levels = month_levels),
    ScanFY        = as.integer(ifelse(ScanMonth > 9,
                                      year(MonthStart) + 1,
                                      year(MonthStart)))
  )

# SubEffort inherits Deploy, Retrieve, UnitType, ScanTimeHrs, DeployedHrs,
# Issue, and Comments directly from ReliableEffort — no second join needed.
SubEffortDetail <- SubEffort |>
  left_join(LocationLabels, by = "Location")

HeatmapData <- SubEffort |>
  left_join(LocationLabels, by = "Location") |>
  group_by(Location, LocationLabel, ScanFY, ScanMonth, ScanMonthName) |>
  summarise(TotalHrs = sum(SubEffortHrs, na.rm = TRUE), .groups = "drop") |>
  # Encode a click key: "FY_MonthNum"
  mutate(CellKey = paste(ScanFY, ScanMonth, sep = "_"))

location_choices <- setNames(LocationLabels$Location, LocationLabels$LocationLabel)

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- page_sidebar(
  title = "PIT Scanner Effort Explorer",
  theme = bs_theme(bootswatch = "cosmo", version = 5),

  sidebar = sidebar(
    width = 260,
    selectInput(
      "location",
      "Location",
      choices  = location_choices,
      selected = location_choices[[1]]
    ),
    hr(),
    p(class = "text-muted small",
      "Click any cell in the heatmap to view the individual deployment records",
      "that contributed to that month's scan hours.")
  ),

  layout_column_wrap(
    width   = 1/2,
    fill    = FALSE,
    heights_equal = "row",
    uiOutput("vbox_deployments"),
    uiOutput("vbox_hours")
  ),

  layout_columns(
    col_widths = 12,
    fill = FALSE,

    card(
      full_screen = TRUE,
      card_header("Reliable Scan Hours by Month and Fiscal Year"),
      plotlyOutput("heatmap", height = "380px")
    ),

    card(
      full_screen = TRUE,
      card_header(
        uiOutput("detail_header")
      ),
      DTOutput("detail_table")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Location-level summary for value boxes
  loc_summary <- reactive({
    ReliableEffort |>
      filter(Location == input$location) |>
      summarise(
        Deployments = n(),
        TotalScanHrs = sum(ScanTimeHrs, na.rm = TRUE)
      )
  })

  output$vbox_deployments <- renderUI({
    s <- loc_summary()
    value_box(
      title = "Reliable Deployments",
      value = format(s$Deployments, big.mark = ","),
      theme = "primary"
    )
  })

  output$vbox_hours <- renderUI({
    s <- loc_summary()
    value_box(
      title = "Total Scan Hours",
      value = format(round(s$TotalScanHrs), big.mark = ","),
      theme = "info"
    )
  })

  # Heatmap data for the selected location
  loc_heatmap <- reactive({
    HeatmapData |> filter(Location == input$location)
  })

  output$heatmap <- renderPlotly({
    df <- loc_heatmap()

    p <- ggplot(df,
                aes(x     = ScanMonthName,
                    y     = factor(ScanFY),
                    fill  = TotalHrs,
                    key   = CellKey,
                    text  = paste0("<b>", ScanMonthName, " FY", ScanFY, "</b>",
                                   "<br>Scan hours: ", round(TotalHrs, 1)))) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_viridis_c(name = "Scan\nHours", option = "plasma",
                           direction = -1, na.value = "grey92") +
      labs(x = NULL, y = "Fiscal Year") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid  = element_blank())

    ggplotly(p, tooltip = "text", source = "effort_heatmap") |>
      layout(clickmode = "event+select") |>
      config(displayModeBar = FALSE)
  })

  # Parse the clicked cell key ("FY_MonthNum") → list(fy, month_num)
  clicked <- reactive({
    ev <- event_data("plotly_click", source = "effort_heatmap")
    req(!is.null(ev), !is.null(ev$key))
    parts <- strsplit(ev$key, "_")[[1]]
    list(fy = as.integer(parts[1]), month_num = as.integer(parts[2]))
  })

  # Detail table: deployment records for the clicked cell
  detail_data <- reactive({
    cl  <- clicked()
    loc <- input$location

    SubEffortDetail |>
      filter(Location == loc,
             ScanFY    == cl$fy,
             ScanMonth == cl$month_num) |>
      arrange(Deploy) |>
      transmute(
        EID,
        `Unit Type`         = UnitType,
        `Deploy`            = format(as.POSIXct(Deploy),   "%Y-%m-%d %H:%M"),
        `Retrieve`          = format(as.POSIXct(Retrieve), "%Y-%m-%d %H:%M"),
        `Sub-Effort Hrs`    = round(SubEffortHrs,  1),
        `Full Deploy Hrs`   = round(DeployedHrs,   1),
        `Scan Time Hrs`     = round(ScanTimeHrs,   1),
        Issue,
        Comments
      )
  })

  output$detail_header <- renderUI({
    cl <- tryCatch(clicked(), error = function(e) NULL)
    if (is.null(cl)) {
      "Deployment Records — click a cell above"
    } else {
      loc_label <- names(location_choices)[location_choices == input$location]
      month_name <- month.abb[cl$month_num]
      paste0("Deployment Records — ", loc_label, "  |  ",
             month_name, " FY", cl$fy)
    }
  })

  output$detail_table <- renderDT({
    cl <- tryCatch(clicked(), error = function(e) NULL)
    if (is.null(cl)) {
      # Return an empty placeholder table
      data.frame(Message = "Click a cell in the heatmap to load records.")
    } else {
      detail_data()
    }
  },
  rownames  = FALSE,
  options   = list(
    pageLength  = 15,
    scrollX     = TRUE,
    dom         = "tip",
    autoWidth   = TRUE,
    columnDefs  = list(list(className = "dt-right",
                            targets   = c(4, 5, 6)))
  ))
}

shinyApp(ui, server)
