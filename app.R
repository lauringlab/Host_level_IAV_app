#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# SETUP
# 
# ggplot()+geom_point(data=f,aes(x=day,y=freq,color=sample_class)) +
# geom_vline(xintercept =
#              unique(f$transmission-f$Donor_onset),linetype=2) + 
#   geom_segment(data=small.min,aes(x=Donor_home_collect-Donor_onset,
#                                   xend=Donor_clinic_collect-Donor_onset,
#                                   y=Donor_home_freq,yend=Donor_clinic_freq),
#                color=cbPalette[5])+
#   geom_segment(data=small.min,aes(x=Recipient_home_collect-Donor_onset,
#                                   xend=Recipient_clinic_collect-Donor_onset,
#                                   y=Recipient_home_freq,yend=Recipient_clinic_freq),
#                color=cbPalette[2])+
#   
#   geom_segment(data=small.min,aes(x=Donor_home_collect-Donor_onset,
#                                   xend=Recipient_home_collect-Donor_onset,
#                                   y=Donor_home_freq,yend=Recipient_home_freq),
#                color=cbPalette[3],linetype=2,alpha=0.5)
# Big function from summary_stats


library("shiny")
library("ggplot2")
require("magrittr")
require(dplyr)
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

##### Transmission
require(readr)
long<-read_csv("./freq.long.csv")
long<- long %>%  mutate(day = collection-Donor_onset,
                        sample_class=gsub(pattern = "(.*)_.*",x=sample,replacement = "\\1",perl = T),
                        pair = paste0(Donor_ENROLLID,"-",Recipient_ENROLLID)) 
wide<-read_csv("./freq.wide.csv")
wide<-wide %>% rowwise() %>% mutate(pair = paste0(Donor_ENROLLID,"-",Recipient_ENROLLID)) 
Houses = setNames(unique(wide$HOUSE_ID),unique(wide$HOUSE_ID))

meta<-read_csv("../Host_level_IAV_evolution/data/processed/secondary/meta_snv_qual.csv")

###################################################################################
# UI
ui <- navbarPage("My Application",
             tabPanel("iSNV overview",
                      pageWithSidebar( 
                        headerPanel("Diversity in Samples"), 
                        
                        sidebarPanel( 
                          HTML("Add explanation of data here"), 
                          width = 3 
                        ), 
                        
                        mainPanel( 
                          
                          # this is an extra div used ONLY to create positioned 
                          # ancestor for tooltip 
                          # we don't change its position 
                          div( 
                            style = "position:relative", 
                            plotOutput("genomePlot",  
                                       hover = hoverOpts("plot_hover", delay = 100,
                                                         delayType = "debounce"), 
                                       dblclick = "plot1_dblclick", 
                                       brush = brushOpts( 
                                         id = "plot1_brush", 
                                         resetOnNew = TRUE)), 
                            uiOutput("hover_info") 
                          ), 
                          width = 7))), 
                          
             tabPanel("Transmission",
                      pageWithSidebar( 
                        headerPanel("Transmission"), 
                        
                        sidebarPanel( 
                            selectInput("select", label = h3("Select House Id"), 
                                        choices = Houses, 
                                        selected = Houses[1]),

                            uiOutput("Pairings")
                            ), 

                        mainPanel( 
                          
                          # this is an extra div used ONLY to create positioned 
                          # ancestor for tooltip 
                          # we don't change its position 
                          div( 
                            style = "position:relative", 
                            plotOutput("transPlot",  
                                       hover = hoverOpts("transPlot_hover", delay = 100,
                                                         delayType = "debounce"), 
                                       dblclick = "trans_dblclick", 
                                       click = "transPlot_click",
                                       brush = brushOpts( 
                                         id = "transPlot_brush", 
                                         resetOnNew = TRUE)), 
                            uiOutput("transPlot_hover_info") 
                          ), 
                          width = 7)
                        )
                      )
             
  )
  

server <- function(input, output) {
  ranges <- reactiveValues(x = NULL, y = NULL)
  ranges_trans <- reactiveValues(x = NULL, y = NULL)
  output$genomePlot <- renderPlot({
    ggplot(figure1,aes(x=pos,y=chr))+
      geom_point(aes(color=class_factor),shape=108,size=10)+
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
                    "<b> Frequency: </b>", round(point$freq.var,2), "<br/>")))
    )
  })
######## transmission tab ##########

  output$Pairings <- renderUI({
    options <- long %>% filter(HOUSE_ID==input$select) %>% 
      select(ends_with("ENROLLID")) %>% 
      mutate(opts = paste0(Donor_ENROLLID,"-",Recipient_ENROLLID)) %>% 
      distinct()
    radioButtons("pairs", "Choose Pair", options$opts)
  })
  
  output$transPlot <- renderPlot({
   long_data<-filter(long,HOUSE_ID==input$select,pair==input$pairs)
   wide_data<-filter(wide,HOUSE_ID==input$select,pair==input$pairs)
    Donor_SPECID<-HIVEr::get_close(meta,date = unique(wide_data$transmission),
                                    enrollid = unique(wide_data$Donor_ENROLLID),
                                    case="donor")
    Recipient_SPECID<-HIVEr::get_close(meta,date = unique(wide_data$transmission),
                                 enrollid = unique(wide_data$Recipient_ENROLLID),
                                 case="recipient")
    select_points<- nearPoints(long_data, input$transPlot_click, allRows = TRUE) # selecting points
    print(select_points$selected_)
    if(any(select_points$selected_)){ # if any are selected dim the others
      selected<-filter(select_points,selected_==T)
      long_data$alpha=0.2
      long_data$alpha[which(long_data$mutation %in% selected$mutation)]<-1
    } else{
      long_data$alpha=1
    }
   if(!(is.na(unique(wide_data$Donor_clinic))) & unique(wide_data$Donor_clinic)==Donor_SPECID){
     Donor_column = "Donor_clinic_freq"
     Donor_time = "Donor_clinic_collect"
   }else if(!(is.na(unique(wide_data$Donor_home))) & unique(wide_data$Donor_home)==Donor_SPECID){
     Donor_column = "Donor_home_freq"
     Donor_time = "Donor_home_collect"
   }
    if(!(is.na(unique(wide_data$Recipient_clinic))) & unique(wide_data$Recipient_clinic)==Recipient_SPECID){
      Recipient_column = "Recipient_clinic_freq"
      Recipient_time = "Recipient_clinic_collect"
    }else if(!(is.na(unique(wide_data$Recipient_home))) & unique(wide_data$Recipient_home)==Recipient_SPECID){
      Recipient_column = "Recipient_home_freq"
      Recipient_time = "Recipient_home_collect"
    }
  ggplot()+geom_point(data=long_data,aes(x=day,y=freq,color=as.factor(sample_class),alpha=alpha)) + 
    scale_color_manual(values = cbPalette[c(5,2)],labels = c("Donor","Recipient"),name="")+ # points (when were the mutations found)
     geom_vline(xintercept =
                  unique(long_data$transmission-long_data$Donor_onset),linetype=2) + # Estimated transmission
     geom_segment(data=wide_data,aes(x=Donor_home_collect-Donor_onset, # lines for donor
                                     xend=Donor_clinic_collect-Donor_onset,
                                     y=Donor_home_freq,yend=Donor_clinic_freq),
                  color=cbPalette[5])+
     geom_segment(data=wide_data,aes(x=Recipient_home_collect-Donor_onset,  # Lines for recipient
                                     xend=Recipient_clinic_collect-Donor_onset,
                                     y=Recipient_home_freq,yend=Recipient_clinic_freq),
                  color=cbPalette[2])+
      geom_segment(data=wide_data,aes(x=get(Donor_time)-Donor_onset,
                                   xend=get(Recipient_time)-Donor_onset,
                                    y=get(Donor_column),
                                    yend=get(Recipient_column)),
                 color=cbPalette[3],linetype=2,alpha=0.5)+
     scale_x_continuous(limits = c(-0.05,max(long_data$day)+0.5))+
     scale_y_continuous(limits = c(min(long_data$freq)-0.02,max(long_data$freq)+0.05))+
     xlab("Days post Donor symptom onset")+ylab("Frequency")+
    scale_alpha_continuous(guide=F)+
      coord_cartesian(xlim = ranges_trans$x, ylim = ranges_trans$y, expand = FALSE)
  })
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$trans_dblclick, {
    brush <- input$transPlot_brush
    if (!is.null(brush)) {
      ranges_trans$x <- c(brush$xmin, brush$xmax)
      ranges_trans$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_trans$x <- NULL
      ranges_trans$y <- NULL
    }
  })
  
  
  output$transPlot_hover_info <- renderUI({
    long_data<-filter(long,HOUSE_ID==input$select,pair==input$pairs)
    hover <- input$transPlot_hover
    point <- nearPoints(long_data, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
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
      p(HTML(paste0("<b> Mutation: </b>", point$mutation, "<br/>",
                    "<b> Frequency: </b>", round(point$freq,2), "<br/>",
                    "<b> Sample: </b>", gsub(pattern = "_",x = point$sample,replacement = ":"), "<br/>"))
    ))
  })
  
  
  
}






runApp(list(ui = ui, server = server))


