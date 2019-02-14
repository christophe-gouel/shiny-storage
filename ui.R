# ui.R

fluidPage(

  titlePanel("Rational expectations storage model"),
  withMathJax(),
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      h2("Parameters choice"),
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "elastD",
                  label = "Demand elasticity, \\(\\alpha^D\\)",
                  min = -0.5,
                  max = -0.05,
                  value = -0.1,
                  round = -2),
      sliderInput(inputId = "k",
                  label = "Storage cost (% of steady-state price, \\(\\bar{P}\\)), \\(k\\)",
                  min = 0,
                  max = 30,
                  value = 5,
                  round = -1),
      sliderInput(inputId = "SDe",
                  label = "Coefficient of variation of the supply shock (%), \\(\\sigma_{e}\\)",
                  min = 0.1,
                  max = 10,
                  value = 2,
                  round = -1),
      submitButton("Update View", icon("refresh")),
      h2("Model definition"),
      helpText("$$\\text{Consumption: } D(P_t) = \\bar{D}\\left[1+\\alpha^{D}\\left(\\frac{P_t-\\bar{P}}{\\bar{P}}\\right)\\right]$$"),
      helpText("$$\\text{Availability: } A_t = S_{t-1}+\\bar{D}\\left(1+\\sigma_{e} e_t\\right)$$"),
      helpText("$$\\text{Storage: } \\beta E_t P_{t+1}-P_t-k\\bar{P}\\le 0, =0 \\text{ if } S_t >0$$"),
      helpText("$$\\text{Market clearing: } A_t = D\\left(P_t\\right)+S_t$$"),
      h3("Fixed parameters"),
      helpText("$$\\bar{D}=\\bar{P}=1$$"),
      helpText("$$\\beta=1/\\left(1+r\\right)=0.98$$"),
      helpText("$$e\\sim N\\left(0,1\\right)$$")
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "demandFun"),
      plotOutput(outputId = "PriceDistrib"),
      h3(textOutput("statHeader")),
      tableOutput("statTable")
    )
  )
)
