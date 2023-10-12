# server.R

function(input, output) {
  model <- reactive({
    model <- list()
    model$params <- list(
      pbar = input$pbar, dbar = input$dbar, k = input$k / 100,
      δ = 0, r = 0.02, elastD = input$elastD, SDe = input$SDe / 100
    )
    model$S <- matrix(
      seq(
        from = 0,
        to = .5 * model$params$dbar,
        length = 200
      ),
      ncol = 1
    )
    model$A <- model$S / 0.5 + model$params$dbar
    model$shocks <- list(
      e = matrix(gherm$nodes, nrow = 1),
      w = matrix(gherm$weights / sum(gherm$weights), ncol = 1)
    )
    model$P <- model$params$pbar * (model$A / model$params$dbar)^model$params$elastD
    model <- SolveStorage(model)
    model$sim <- SimulateStorage(model, 1)
    return(model)
  })


  output$demandFun <- renderPlot({
    m <- model()
    validate(need(m$SolveStat$exitflag == 1, "Failure to solve the model"))
    A <- seq(
      from = m$params$dbar * exp(-4 * m$params$SDe),
      to = max(m$A),
      length = 1000
    )
    dt <- tibble(
      A = A,
      Pinv = m$params$pbar * (A / m$params$dbar)^(1 / m$params$elastD),
      P = m$PriceFunction(A)
    )

    dt %>%
      gather(key = Variable, value = value, -A) %>%
      ggplot(aes(x = A, y = value, col = Variable)) +
      geom_line() +
      scale_color_manual(
        values = c("red", "black"),
        name = "",
        breaks = c("Pinv", "P"),
        labels = c(
          "Inverse demand function",
          "Price including demand for storage"
        )
      ) +
      geom_point(aes(x = m$params$dbar, y = m$params$pbar, col = "Pinv")) +
      xlab("Availability") +
      ylab("Price") +
      ggtitle("Inverse demand function") +
      theme_light() +
      theme(
        legend.position = "top",
        legend.justification = c(1, 1)
      )
  })

  output$PriceDistrib <- renderPlot({
    m <- model()
    validate(need(m$SolveStat$exitflag == 1, "Failure to solve the model"))
    sim <- m$sim
    ggplot(data = tibble(P = c(sim$P)), aes(x = P)) +
      geom_histogram(aes(y = ..density.., color = "P"), colour = "red", fill = "white") +
      geom_density(alpha = .2, fill = "#FF6666", col = "red") +
      geom_line(aes(y = P, x = d),
        data = tibble(
          d = seq(min(c(sim$P)), max(c(sim$P)), len = 100),
          P = dnorm(d,
            mean = m$params$pbar,
            sd = -m$params$pbar * m$params$SDe / m$params$elastD
          )
        )
      ) +
      xlab("Price") +
      ylab("Density") +
      geom_vline(xintercept = mean(c(sim$P)), col = "blue") +
      ggtitle("Price distribution (in black price distribution without stocks, in blue mean price)") +
      theme_light()
  })

  output$statHeader <- renderText("Statistics on the asymptotic distribution")

  output$statTable <- renderTable({
    m <- model()
    validate(need(m$SolveStat$exitflag == 1, "Failure to solve the model"))
    P_acf <-
      m$sim$P %>%
      apply(1, function(x) acf(x, plot = FALSE)[1:2]$acf) %>%
      rowMeans()
    Psim_no <- m$invdemand(m$params$dbar * rlnorm(1E6, sdlog = m$params$SDe))
    A_acf <-
      m$sim$A %>%
      apply(1, function(x) acf(x, plot = FALSE)[1:2]$acf) %>%
      rowMeans()
    tibble(P = c(m$sim$P), A = c(m$sim$A)) %>%
      gather(key = Variable) %>%
      group_by(Variable) %>%
      summarize(
        Mean = mean(value),
        Median = median(value),
        CV = sd(value) / mean(value),
        Skewness = skewness(value),
        `Excess-kurtosis` = kurtosis(value) - 3,
        Q05 = quantile(value, 0.05),
        Q95 = quantile(value, 0.95)
      ) %>%
      mutate(
        AC1 = c(A_acf[1], P_acf[1]),
        AC2 = c(A_acf[2], P_acf[2])
      ) %>%
      add_row(
        Variable = "A (without stocks)",
        Mean = m$params$dbar * exp(m$params$SDe^2 / 2),
        Median = m$params$dbar,
        CV = m$params$dbar * sqrt((exp(m$params$SDe^2) - 1) * exp(m$params$SDe^2)) / Mean,
        Skewness = (exp(m$params$SDe^2) + 2) * sqrt(exp(m$params$SDe^2) - 1),
        `Excess-kurtosis` = exp(4 * m$params$SDe^2) + 2 * exp(3 * m$params$SDe^2) + 3 * exp(2 * m$params$SDe^2) - 6,
        Q05 = m$params$dbar * qlnorm(0.05, sdlog = m$params$SDe),
        Q95 = m$params$dbar * qlnorm(0.95, sdlog = m$params$SDe),
        AC1 = 0, AC2 = 0
      ) %>%
      add_row(
        Variable = "P (without stocks)",
        Mean = mean(Psim_no),
        Median = median(Psim_no),
        CV = sd(Psim_no) / Mean,
        Skewness = skewness(Psim_no),
        `Excess-kurtosis` = kurtosis(Psim_no) - 3,
        Q05 = quantile(Psim_no, 0.05),
        Q95 = quantile(Psim_no, 0.95),
        AC1 = 0, AC2 = 0
      )
  })
}
