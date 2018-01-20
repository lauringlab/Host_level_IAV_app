#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library("shiny")
library("ggplot2")
require(wesanderson)
cbPalette<-wes_palette("Zissou")
theme_set(new = theme_classic()+ theme(
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour ='black',size=0.5,linetype='solid'),
  text=element_text(family="Arial",size = 18))) # to make nice plots
# get the start and stop of each OR for each segment (2014-2015 used as reference)

figure1<-read.csv("../Host_level_IAV_evolution/results/Figures/data/Figure1C_D.csv",
                  stringsAsFactors = F)
# Set the segments as factors with PB2 on top
figure1$chr<-factor(figure1$chr,
                    levels = rev(c("PB2","PB1","PA","HA","NP","NR","M","NS"))) 

chrs<-read.csv("../Host_level_IAV_evolution/data/reference/segs.csv",stringsAsFactors = T) 
chrs$chr<-factor(chrs$chr,levels=levels(figure1$chr)) # set factors on the is meta data
ui <- pageWithSidebar(
  headerPanel("Tooltips in ggplot2 + shiny"),
  
  sidebarPanel(
    HTML("Add explanation of data here"),
    width = 3
  ),
  
  mainPanel(
    
    # this is an extra div used ONLY to create positioned ancestor for tooltip
    # we don't change its position
    div(
      style = "position:relative",
      plotOutput("genomePlot", 
                 hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce"),
                 dblclick = "plot1_dblclick",
                 brush = brushOpts(
                   id = "plot1_brush",
                   resetOnNew = TRUE)),
      uiOutput("hover_info")
    ),
    width = 7
  )
)

server <- function(input, output) {
  ranges <- reactiveValues(x = NULL, y = NULL)

  output$genomePlot <- renderPlot({
    ggplot(figure1,aes(x=pos,y=chr))+
      geom_point(aes(color=class_factor),shape=108,size=8)+
      geom_segment(data=chrs,aes(x = start, y = chr, xend = stop, yend = chr))+
      ylab("")+
      xlab("")+
      scale_color_manual(name="",values=cbPalette[c(1,4)])+
      theme(axis.ticks =element_blank(),
            axis.line.x = element_blank(), axis.line.y=element_blank())+
      scale_x_continuous(breaks=c())+
      theme(legend.position = "none")+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    

  })
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  
  output$hover_info <- renderUI({
    hover <- input$plot_hover
    point <- nearPoints(figure1, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px;")
    
    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b> Isolate: </b>", point$SPECID, "<br/>",
                    "<b> Position: </b>", point$pos, "<br/>",
                    "<b> Frequency: </b>", round(point$freq.var,2), "<br/>",
                    "<b> Distance from left: </b>", left_px, "<b>, from top: </b>", top_px)))
    )
  })
}

runApp(list(ui = ui, server = server))


