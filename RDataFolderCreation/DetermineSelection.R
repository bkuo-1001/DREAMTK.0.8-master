DeterminingSelection = read.csv("C:/ProgramData/MySQL/MySQL Server 5.6/Determiningselection_data.csv",header = TRUE, sep = ";", na.strings = "", stringsAsFactors = FALSE);



select_labels <- sapply(DeterminingSelection$selected, function(x){
  if(x == "1"){
    return ("Selected");
  } else {
    return ("Not Selected");
  }
})
library(ggplot2);
library(plotly);
p<- plot_ly(DeterminingSelection, x = ~aic, y = ~ac50, alpha = 0.5,
                              colors = rainbow(10, s = 3/4, v=3/4),
                              type = "scatter", mode = "markers",
                              symbol = select_labels, symbols = c("circle", "x"),
                              hoverinfo = "text",
                              
                              text = ~paste( "</br> Casn: ", casn,
                                             "</br> Assay endpoint: ", assay_component_endpoint_name,
                                             "</br> Selected: ", selected,
                                             "</br> Ac50 (uM): ", signif(ac50, digits=5),
                                             "</br> AIC: ", signif(aic, digits=5))  ) %>% 
  
  layout(title = "AIC VS Ac50",
         xaxis = list(title = "AIC",
                      titlefont = list(size = 14)),
         yaxis = list(title = "Ac50",
                      titlefont = list(size = 14)),
         showlegend = TRUE,
         legend = list(orientation = "h"))

print(p)
p2<- plot_ly(DeterminingSelection, x = ~prob, y = ~ac50, alpha = 0.5,
            colors = rainbow(10, s = 3/4, v=3/4),
            type = "scatter", mode = "markers",
            symbol = select_labels, symbols = c("circle", "x"),
            hoverinfo = "text",
            
            text = ~paste( "</br> Casn: ", casn,
                           "</br> Assay endpoint: ", assay_component_endpoint_name,
                           "</br> Selected: ", selected,
                           "</br> Ac50 (uM): ", signif(ac50, digits=5),
                           "</br> PROB: ", signif(prob, digits=5))  ) %>% 
  
  layout(title = "PROB VS Ac50",
         xaxis = list(title = "PROB",
                      titlefont = list(size = 14)),
         yaxis = list(title = "Ac50",
                      titlefont = list(size = 14)),
         showlegend = TRUE,
         legend = list(orientation = "h"))
print(p2)
p3<- plot_ly(DeterminingSelection, x = ~prob, y = ~aic, alpha = 0.5,
             colors = rainbow(10, s = 3/4, v=3/4),
             type = "scatter", mode = "markers",
             symbol = select_labels, symbols = c("circle", "x"),
             hoverinfo = "text",
             
             text = ~paste( "</br> Casn: ", casn,
                            "</br> Assay endpoint: ", assay_component_endpoint_name,
                            "</br> Selected: ", selected,
                            "</br> Ac50 (uM): ", signif(ac50, digits=5),
                            "</br> AIC: ", signif(aic, digits=5))  ) %>% 
  
  layout(title = "AIC VS PROV",
         xaxis = list(title = "PROB",
                      titlefont = list(size = 14)),
         yaxis = list(title = "AIC",
                      titlefont = list(size = 14)),
         showlegend = TRUE,
         legend = list(orientation = "h"))
print(p3)