# Change History

## Update: 02/12/2026

Selector tab — Positive and negative lasso selection:                                                                              
                                                                                
  - First lasso (no shift) → sets ingroup, highlighted with muted red shape + "Ingroup" label                                        
  - Shift + lasso → sets outgroup, highlighted with muted blue shape + "Outgroup" label; multiple shift-lassos accumulate            
  - "Find markers" button → compares ingroup vs outgroup if both defined, otherwise ingroup vs all other cells
  - New non-shift lasso → resets both groups
  - Double-click → clears everything
  
  ui.R (2 edits):                                                                                                                    
  - Line 758: Updated tab description to: "Lasso cells on UMAP to discover markers (top 3000 HVGs). Hold Shift and draw a second     
  lasso to define an explicit outgroup."
  - Line 762: Updated hint text to: "Select cells with lasso or box tool. Hold Shift and draw a second selection to define an
  outgroup."

  server.R (4 edits):

  1. Lines 1167-1174 — Replaced selected_cells_rv (which used event_data("plotly_selected")) with two reactives driven by JS-set
  Shiny inputs:
    - selected_cells_rv reads input$sel_ingroup_keys
    - outgroup_cells_rv reads input$sel_outgroup_keys
  2. Lines 1262-1416 — Rewrote renderPlotly for sel_umap:
    - Added dragmode = "lasso" to layout
    - Added htmlwidgets::onRender() JavaScript that:
        - Tracks Shift key state via keydown/keyup listeners
      - On plotly_selected: first lasso → ingroup; Shift+lasso → outgroup (= all newly selected keys minus ingroup)
      - Draws colored shape overlays: muted red (rgba(205,92,92,0.15)) for ingroup, muted blue (rgba(100,149,237,0.15)) for outgroup
      - Adds bold text annotations ("Ingroup" / "Outgroup") positioned at the top of each shape
      - Sends cell IDs to Shiny via Shiny.setInputValue
      - Handles plotly_deselect to clear everything
  3. Lines 1419-1429 — Updated sel_ncells to show:
    - "Ingroup: X cells | Outgroup: Y cells" when outgroup is defined
    - "X cells (vs all other cells)" when only ingroup is selected
  4. Lines 1438-1447 — Modified do_marker handler:
    - If outgroup_cells_rv() has cells → uses them as group2
    - Otherwise falls back to all-other-cells (original behavior)
    
  - Added lastDrawnPath variable and a plotly_selecting listener that captures the lasso/box path during each drag — this fires
  reliably for every draw operation regardless of shift state
  - plotly_selected now uses lastDrawnPath as the primary path source, falling back to the event's own lassoPoints/range only if the
  selecting event didn't fire
  - Renamed getPathInfo → buildPathInfo (now takes explicit args instead of the whole event object) for cleaner reuse between
  plotly_selecting and plotly_selected
  - lastDrawnPath is cleared after use in plotly_selected and on deselect
  
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

