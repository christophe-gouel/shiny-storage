# server.R

function(input, output){

  model <- reactive({
    model <- list()
    model$params <- list(pbar = input$pbar, dbar = input$dbar, k = input$k / 100, delta = 0,
                         r = 0.02, elastD = input$elastD, SDe = input$SDe / 100)
    model$s <- matrix(seq(model$params$dbar - 4 * model$params$SDe,
                          1.5 * model$params$dbar,
                          length = 200), ncol = 1)
    model$shocks <- list(e = matrix(gherm$nodes, nrow = 1),
                         w = matrix(gherm$weights/sum(gherm$weights), ncol = 1))
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
      geom_point(aes(x = m$params$dbar, y = m$params$pbar, col = "Pinv")) +
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
                CV = sd(value) / mean(value),
                Skewness = skewness(value),
                `Excess-kurtosis` = kurtosis(value)-3,
                Q05 = quantile(value, 0.05),
                Q95 = quantile(value, 0.95)) %>%
      mutate(AC1 = c(A_acf[1], P_acf[1]),
             AC2 = c(A_acf[2], P_acf[2])) %>%
      add_row(Variable = "A (without stocks)",
              Mean = m$params$dbar,
              Median = m$params$dbar,
              CV = m$params$dbar*m$params$SDe / m$params$dbar,
              Skewness = 0, `Excess-kurtosis` = 0,
              Q05 = qnorm(0.05, mean = m$params$dbar, sd = m$params$dbar * m$params$SDe),
              Q95 = qnorm(0.95, mean = m$params$dbar, sd = m$params$dbar * m$params$SDe),
              AC1 = 0, AC2 = 0) %>%
      add_row(Variable = "P (without stocks)", Mean = m$params$pbar, Median = m$params$pbar,
              CV = -m$params$pbar*m$params$SDe/m$params$elastD / m$params$pbar,
              Skewness = 0, `Excess-kurtosis` = 0,
              Q05 = qnorm(0.05, mean = m$params$pbar, sd = -m$params$pbar * m$params$SDe / m$params$elastD),
              Q95 = qnorm(0.95, mean = m$params$pbar, sd = -m$params$pbar * m$params$SDe / m$params$elastD),
              AC1 = 0, AC2 = 0)

  })
}

