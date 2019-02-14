# server.R

function(input, output){

  model <- reactive({
    model <- list()
    model$params <- list(pbar = 1, dbar = 1, k = input$k / 100, delta = 0,
                         r = 0.02, elastD = input$elastD, SDe = input$SDe / 100)
    model$s <- seq(0.8, 1.5, length = 200)
    model$shocks <- list(e = gherm$nodes,
                         w = gherm$weights/sum(gherm$weights))
    model$P <- model$params$pbar*(1+(model$s-model$params$dbar)/(model$params$elastD*model$params$dbar))
    model <- SolveStorage(model)
    model$sim <- SimulateStorage(model, 1)
    return(model)
  })


  output$demandFun <- renderPlot({
    m <- model()
    validate(need(m$SolveStat$exitflag == 1, "Failure to solve the model"))
    dt <- tibble(A = m$s,
                 Pinv = m$params$pbar*(1+(A-m$params$dbar)/(m$params$elastD * m$params$dbar)),
                 P = m$P)

    dt %>%
      gather(key = Variable, value = value, -A) %>%
      ggplot(aes(x = A, y = value, col = Variable)) +
      geom_line() +
      scale_color_manual(values = c("red", "black"),
                         name = "",
                         breaks = c("Pinv", "P"),
                         labels = c("Inverse demand function",
                                    "Price including demand for storage")) +
      geom_point(aes(x = 1, y = 1, col = "Pinv")) +
      xlim(min(c(m$sim$A))*0.9, max(c(m$sim$A)*1.1)) +
      ylim(min(c(m$sim$P))*0.9, max(c(m$sim$P)*1.1)) +
      xlab("Availability") + ylab("Price") +
      ggtitle("Inverse demand function") +
      theme_light() +
      theme(legend.position = "top",
            legend.justification = c(1, 1))
  })

  output$PriceDistrib <- renderPlot({
    m <- model()
    validate(need(m$SolveStat$exitflag == 1, "Failure to solve the model"))
    sim <- m$sim
    ggplot(data = tibble(P = c(sim$P)), aes(x = P)) +
      geom_histogram(aes(y = ..density.., color = "P"), colour = "red", fill = "white") +
      geom_density(alpha = .2, fill = "#FF6666", col = "red") +
      geom_line(aes(y = P, x = d),
                data = tibble(d = seq(min(c(sim$P)), max(c(sim$P)), len = 100),
                              P = dnorm(d, mean = m$params$pbar,
                                        sd = -m$params$pbar*m$params$SDe/m$params$elastD))) +
      xlab("Price") + ylab("Density") +
      geom_vline(xintercept = mean(c(sim$P)), col = "blue") +
      ggtitle("Price distribution (in black price distribution without stocks, in blue mean price)") +
      theme_light()
  })

  output$statHeader = renderText("Statistics on the asymptotic distribution")

  output$statTable <- renderTable({
    m <- model()
    validate(need(m$SolveStat$exitflag == 1, "Failure to solve the model"))
    P_acf <-
      m$sim$P %>%
      apply(1, function(x) acf(x, plot = FALSE)[1:2]$acf) %>%
      rowMeans()
    A_acf <-
      m$sim$A %>%
      apply(1, function(x) acf(x, plot = FALSE)[1:2]$acf) %>%
      rowMeans()
    tibble(P = c(m$sim$P), A = c(m$sim$A)) %>%
      gather(key = Variable) %>%
      group_by(Variable) %>%
      summarize(Mean = mean(value),
                Median = median(value),
                Std = sd(value),
                Skewness = skewness(value),
                `Excess-kurtosis` = kurtosis(value)-3,
                Q05 = quantile(value, 0.05),
                Q95 = quantile(value, 0.95)) %>%
      mutate(AC1 = c(A_acf[1], P_acf[1]),
             AC2 = c(A_acf[2], P_acf[2])) %>%
      add_row(Variable = "P (without stocks)", Mean = 1, Median = 1,
              Std = -m$params$pbar*m$params$SDe/m$params$elastD,
              Skewness = 0, `Excess-kurtosis` = 0,
              Q05 = qnorm(0.05, mean = 1, sd = -m$params$pbar*m$params$SDe/m$params$elastD),
              Q95 = qnorm(0.95, mean = 1, sd = -m$params$pbar*m$params$SDe/m$params$elastD),
              AC1 = 0, AC2 = 0)

  })
}
