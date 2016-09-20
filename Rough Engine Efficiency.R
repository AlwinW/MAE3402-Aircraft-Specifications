# NB requires functions from "Helper Calculation Functions.R"

RepeatRows <- function(data, each = 1) {
  data <- data[rep(row.names(data), each = each), 1:length(data)]
  return(as.data.frame(data))
}

UpdateParams <- function(input) {
  input$b <- sqrt(input$AR * input$S)
  input$K <- 1/(pi * input$AR * input$e)
  input$W <- input$m * 9.8065
  input$WS <- input$W / input$S
  input$P0 <- input$P0eng * 2
  return(input)
}

WS <- seq(1500, 2900, by = 10)
varyWS <- RepeatRows(input_initial, length(WS))
varyWS$WS <- WS
varyWS <- mutate(varyWS, m = WS * S / 9.8065) %>%
  UpdateParams(.) %>%
  InpSpecs(., specifications) 

varyWS <- varyWS %>%
  mutate(TOP = W/P0 * WS * 1/(Clclean + Clflaps),
         TOFL = 11.8 * TOP + 0.255 * TOP^2)

TakeOffLength = data.frame(WS = varyWS$WS, 
                           NTO = NA,
                           BFL = NA,
                           IntegralApproach = NA,
                           TOPApproach = varyWS$TOFL)



for (i in 1:nrow(varyWS)) {
  temp = TakeOff(varyWS[i,])
  TakeOffLength$NTO[i] = temp$NTO
  TakeOffLength$BFL[i] = temp$BFL
  TakeOffLength$IntegralApproach[i] = temp$TakeOffDistance
}
   

ggplot(TakeOffLength) + 
  geom_line(aes(x = WS, y = NTO, colour = "NTO")) + 
  geom_line(aes(x = WS, y = BFL, colour = "BFL")) 

ggplot(TakeOffLength) + 
  geom_line(aes(x = WS, y = IntegralApproach, colour = "IntegralApproach")) + 
  geom_line(aes(x = WS, y = TOPApproach, colour = "TOPApproach"))

ggplot(TakeOffLength) +
  geom_line(aes(x = WS, y = TOPApproach / IntegralApproach, colour = "TOP/Int"))

ggplot(TakeOffLength) +
  geom_line(aes(x = WS, y = TOPApproach / IntegralApproach - lag(TOPApproach, 1) / lag(IntegralApproach, 1), colour = "TOP/Int"))
