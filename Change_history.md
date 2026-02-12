# Change History

## Update: 02/11/2026

  ui.R (lines 564-572)

  Added a new optional dropdown sc1c1inp1b ("Cell information (X-axis fill, optional)") that allows users to select a secondary
  categorical variable for double grouping. The default is "(none)" for single grouping behavior.

  server.R - scVioBox function (lines 389-468)

  Modified the function to:
  - Accept a new inp1b parameter for the secondary grouping variable
  - Check if secondary grouping is enabled (useFill)
  - When enabled: add the fill variable as X2, use aes(X, val, fill = X2), and display a legend
  - When disabled: maintain original behavior with aes(X, val, fill = X) and no legend

  server.R - Function calls (lines 1031, 1044, 1054)

  Updated all three calls to scVioBox (renderPlot, PDF download, PNG download) to pass input$sc1c1inp1b.

  How it works

  - Single grouping (default): User selects only "Cell information (X-axis)" - plot shows violins/boxplots grouped and colored by
  that variable (no legend)
  - Double grouping: User also selects "Cell information (X-axis fill, optional)" - plot shows violins/boxplots grouped by the first
  variable on X-axis, but colored/filled by the second variable (with legend)

