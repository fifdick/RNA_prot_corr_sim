library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggnewscale)
library(extraDistr)
# Let's assume perfect correlation between RNA and protein for this "toy" gene
# Let's assume gene x is expressed only in the some and translated in the soma and x % the protein products are transported to synapses in YG and CT but x % in PD
# Let's assume the true RNA expression is median(RNA_cell_norm)
# Protein expression calculated on this assumption of true RNA expression
# RNA protein correlation based on the RNA expression drawn from the distribution and normalized by library size


set.seed(1)


# soma/synapse ratio
# 0.5 means as much soma as synapse
# < 0.5 means more synapses than soma
# > 0.5 means more soma than synapses
# all samples (across groups) get r pulled from the same distribition
#(although later we could shift it a bit for YG)
#lets assume 100% of gene expr is in soma and 0% protein, 100% in synapse and 0% gene expr in synapse
#lets assume perfect RNA - prot correlation

library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("RNA-Prot correlation in neurons"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the number of bins ----
      sliderInput(inputId = "sample_size",
                  label = "Number of samples per group:",
                  min = 10,
                  max = 100,
                  value = 20),
     sliderInput(label = "Bin width in histogram",
                  inputId = "bins",
                  min = 0.01,
                  max = 0.2,
                  value = 0.1),
  # Input: Slider for the number of bins ----
      sliderInput(inputId = "var_libSize",
                  label = "Variance of RNA library size across samples:",
                  min = 0,
                  max = 100000,
                  value = 10000),
  # Input: Slider for the number of bins ----
      sliderInput(inputId = "perc_g_expr",
                  label = "Gene expression (% of libsize)",
                  min = 0.0000001,
                  max = 0.0001,
                  value = 0.0001),
      selectInput(
  		inputId = "sizep",
  		label = "Dispersion gene expression",
  		choices = c("0.001", "0.01", "0.1", "0.3"),
  		selected = 0.01,
  		multiple = FALSE,
  		selectize = TRUE,
  		width = NULL,
  		size = NULL
		),
      sliderInput(inputId = "perc_RNA_soma",
                  label = "Percentage of RNA residing in soma",
                  min = 0,
                  max = 1,
                  value = 0.9),
     sliderInput(inputId = "perc_prot_soma",
                  label = "Percentage of protein residing in soma",
                  min = 0,
                  max = 1,
                  value = 0.1),
    hr(),
      checkboxInput("r_const", label = "Constant soma:synapse ratio", value = FALSE),
      sliderInput(inputId = "r_const_val",
		  label = "Value of soma:synapse ratio (if constant)",
		  min = 0,
		  max = 1,
		  value = 0.5),
    hr()
      ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "corrPlot", width = 500, height = 700),
      # HTML('<img src="img1.png", height="400px"    
       #   style="float:right"/>','<p style="color:black"></p>')
      img(src = "img1.png", width = 500, height = 500)

    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output, session) {
 observe({
        val <- input$perc_RNA_soma
        # Control the value, min, max, and step.
        # Step size is 2 when input value is even; 1 when value is odd.
        updateSliderInput(session, "perc_prot_soma",
			  value = 1-val,
          		  min = 0,
			  max = 1)
      })
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$corrPlot <- renderPlot({

    	n = input$sample_size
	#soma_syn_r <- rbeta(n = n*3, shape1 = 5, shape2 = 5)
	set.seed(1)
	if(input$r_const == FALSE){
		soma_syn_r<- c(rbeta(n = n, shape1 = 8, shape2 = 2), 
		     rbeta(n = n, shape1 = 5, shape2 = 5),
		     rbeta(n = n, shape1 = 2, shape2 = 8))
  	} else {
	        soma_syn_r <- c(rep(input$r_const_val, n*3))			
	}


	RNA_libsize  <- rnorm(n*3, 40000000, input$var_libSize) 

	#The negative binomial distribution, especially in its alternative parameterization described above, can be used as an alternative to the Poisson distribution. It is especially useful for discrete data over an unbounded positive range whose sample variance exceeds the sample mean. In such cases, the observations are overdispersed with respect to a Poisson distribution, for which the mean is equal to the variance. Hence a Poisson distribution is not an appropriate model. Since the negative binomial distribution has one more parameter than the Poisson, the second parameter can be used to adjust the variance independently of the mean. See Cumulants of some discrete probability distributions.
	#In addition to James MacDonald's answer, here's another relevant passage from the same Wikipedia article: "The negative binomial distribution also arises as a continuous mixture of Poisson distributions (i.e. a compound probability distribution) where the mixing distribution of the Poisson rate is a gamma distribution." The other important bit of information to know is that read counts for a sample in theory follow a Binomial(n,p) distribution, where n is the total number of reads and p is the probability of a read mapping to a specific gene. However, the binomial distribution is computationally inconvenient, since it involves computing factorials, and with millions of reads in each sample, the Poisson distribution (with lambda = n*p) is an excellent approximation to the Binomial while being far more mathematically tractable. So the Poisson noise quantifies the counting uncertainty, while the gamma distribution quantifies the variation in gene expression between replicates. The mixture of the two yields the negative binomial. See also section 2 of http://www.statsci.org/smyth/pubs/edgeRChapterPreprint.pdf

  # TODO: think about better parameters of the nb distr
  # or make a slider for them to see how they effect the model
	#RNA_expr  <-  ((1000000/100) * input$perc_g_expr) / RNA_libsize
	RNA_expr_counts_unn <- rnbinom(n = n*3, size = 1/as.numeric(input$sizep), mu = input$perc_g_expr * 40000000) 
	RNA_expr <- RNA_expr_counts_unn / RNA_libsize 

	#define the percentage of RNA expr in the soma
	#expr_soma = input$perc_RNA_soma # prot expr soma = 1-expr_soma

	RNA_tissue  <-  ((soma_syn_r) * RNA_expr * input$perc_RNA_soma) + (1-soma_syn_r) * RNA_expr * (1-input$perc_RNA_soma) 
	Prot_tissue  <- (soma_syn_r * input$perc_prot_soma * RNA_expr) + ((1-input$perc_prot_soma) * RNA_expr * (1-soma_syn_r))

	df <- data.frame(RNA = RNA_tissue,
		 Protein = Prot_tissue,
		 somaTOsynapse = soma_syn_r,
		condition = c(rep("YG", n), rep("HA", n), rep("PD", n)),
		RNA_expr_counts_unn = RNA_expr_counts_unn)

	cols <- c("PD" = "#ef476f", "HA" = "#073b4c", "YG" = "#06d6a0")

	ggplot(df, aes(x = RNA, y = Protein)) +
	geom_smooth(method = "lm", col = "lightgrey") +
	geom_point(aes(col = somaTOsynapse)) +		
	scale_color_gradient(high = "#bc9504", low = "#009c9c") + 
	facet_wrap(~condition) +
	new_scale_color() +
	ggpubr::stat_cor(aes(col = condition)) +
	scale_color_manual(values = cols) +
	theme_light() -> p1
	
	ggplot(df,aes(x = soma_syn_r)) +
		geom_density(aes(col = condition)) +
		geom_histogram(aes(y = ..density.., fill = condition), alpha = .1, binwidth = input$bins) +
		theme_light() +
		scale_fill_manual(values = cols) +
		scale_color_manual(values = cols) + 
		guides(fill = F) -> p2
	ggplot(df, aes(x = RNA_expr_counts_unn)) +
		geom_density() -> p3
	ggplot(df, aes(x = RNA_expr)) +
		geom_density() -> p4
	p1 / (plot_spacer() +  p2) / (p3 + p4)# + patchwork::plot_layout(widths=c(1,0.3))

   
    })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
