source(file.path("utils", "helpers.R"),  local = TRUE)

# Used to dynamically generate selectable organisms based on taxlev
output$sra_order_organisms <- renderUI({
  tables <- pstat.extraction(pstat)
  TAX_TABLE <- tables$TAX
  choices <- unique(TAX_TABLE[[input$sra_taxlev]])
  selectizeInput('sra_order_organisms', label='Reorder Organisms', choices=choices, multiple=TRUE)
})

output$ra_plot <- renderPlotly({

  # Future proofing this plot
  tables <- pstat.extraction(pstat)
  TAX_TABLE <- tables$TAX
  RA_TABLE <- tables$RLA
  SAM_DATA <- tables$SAM

  # Sum by taxon level
  df.ra <- upsample.ra(RA_TABLE, TAX_TABLE, input$sra_taxlev)

  # If grouping is selected
  if (input$group_samples & !is.null(input$gra_select_conditions)) {
    df.ra$covariate <- SAM_DATA[[input$gra_select_conditions]]
    df.ra.melted <- melt(df.ra, id.vars = "covariate")
    df.avg.ra <- aggregate( . ~ variable + covariate , data = df.ra.melted, mean)
    df.avg.ra <- dcast(data = df.avg.ra,formula = covariate~variable)
    rownames(df.avg.ra) <- df.avg.ra$covariate
    df.sam <- df.avg.ra[,"covariate",drop=FALSE]
    colnames(df.sam) <- input$gra_select_conditions
    df.avg.ra$covariate <- NULL
    df.avg.ra <- df.avg.ra[,order(colSums(df.avg.ra))]
    df.ra <- df.avg.ra
  }

  # Reorder by most prominent organisms
  df.ra <- df.ra[,order(colSums(df.ra))]

  # Put selected organisms first
  if (!is.null(input$sra_order_organisms)) {
    organisms.order <- c(setdiff(colnames(df.ra), input$sra_order_organisms), rev(input$sra_order_organisms))
    df.ra <- df.ra[,organisms.order]
  }

  # Order samples by organisms if not by conditons
  if (input$sra_sort_by == "organisms") {
    for (i in 1:ncol(df.ra)) {
      df.ra <- df.ra[order(df.ra[,i]),]
    }
  }

  # If any conditions are selected make a side bar
  if (!is.null(input$sra_select_conditions) || input$group_samples) {
    
    if (!input$group_samples) {
      df.sam <- SAM_DATA[,input$sra_select_conditions,drop=FALSE]
    }

    # Order samples by conitions if not by organisms
    if (input$sra_sort_by == "conditions") {
      for (i in ncol(df.sam):1) {
        df.sam <- df.sam[order(df.sam[[i]]),,drop=FALSE]
      }
      # Reorder stacked barplot
      df.ra <- df.ra[order(match(rownames(df.ra), rownames(df.sam))),,drop=FALSE]
    } else {
      df.sam <- df.sam[order(match(rownames(df.sam), rownames(df.ra))),,drop=FALSE]
    }

    # Retain hover-text information before conditions are factorized
    hover.txt <- c()
    for (i in 1:ncol(df.sam)) {
      hover.txt <- cbind(hover.txt, df.sam[[i]])
    }

    # Plotly | Heatmap
    df.sam[] <- lapply(df.sam, factor)
    m <- data.matrix(df.sam)
    m.row.normalized <- apply(m, 2, function(x)(x-min(x))/(max(x)-min(x)))
    hm <- plot_ly(x = colnames(m), y = rownames(m), z = m.row.normalized, 
                  type = "heatmap",
                  showscale=FALSE,
                  hoverinfo = "x+y+text",
                  text=hover.txt) %>%
           layout(xaxis = list(title = "", tickangle = -45),
                  yaxis = list(showticklabels = FALSE, type = 'category', ticks = ""))
  }

  # Plotly | Stacked Bar Plots
  df.plot <- df.ra
  df.plot$samples <- rownames(df.plot)
  sbp <- plot_ly(df.plot, y = ~samples, x = df.plot[[colnames(df.plot)[1]]], 
                 type = 'bar', 
                 orientation = 'h', 
                 name = substr(colnames(df.plot)[1], 1, 40)) %>%
          layout(font = list(size = 10),
                 yaxis = list(title = '', type = 'category',
                              tickmode = "array",
                              tickvals = rownames(df.plot),
                              showticklabels = FALSE,
                              categoryorder = 'trace'),
                 xaxis = list(title = 'Relative Abundance'),
                 barmode = 'stack',
                 showlegend = input$sra_show_legend)
  for (i in 2:(ncol(df.plot)-1)) {
    sbp <- add_trace(sbp, x=df.plot[[colnames(df.plot)[i]]], name=substr(colnames(df.plot)[i], 1, 40))
  } 

  # Create a multiplot if any conditions are selected
  if (!is.null(input$sra_select_conditions) || input$group_samples) {
    hm.sbp <- subplot(hm, sbp, widths=c(0.1,  0.9))
    hm.sbp$elementId <- NULL # To suppress a shiny warning
    return(hm.sbp)
  } else {
    sbp$elementId <- NULL # To suppress a shiny warning
    return(sbp)
  }
})