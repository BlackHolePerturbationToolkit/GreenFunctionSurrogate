(* ::Package:: *)

Paclet[
  Name -> "GreenFunctionSurrogate",
  Version -> "1.0.0",
  MathematicaVersion -> "10+",
  Creator -> "Barry Wardell",
  Description -> "A set of functions for evaluating a surrogate model for the Green function for the Regge-Wheeler equation.",
  Extensions ->
  {
    { "Kernel",
	    "Context" -> {
        "GreenFunctionSurrogate`",
        "GreenFunctionSurrogate`Evaluate`"
      }
	  },

    { "Documentation",
      Language -> "English", 
      MainPage -> "Guides/GreenFunctionSurrogate",
      Resources -> 
     	{
        "Guides/GreenFunctionSurrogate"
      }
    }
  }
]
