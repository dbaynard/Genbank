module Bio.GenbankParser.Combinators (
    eol
  , whiteSpace
  , notChar
  , noEol
  , noWhiteSpace
  , space
  , spaces
  , spaces1
)	where

import Text.ParserCombinators.Parsec hiding (space, spaces)
{-import Text.ParserCombinators.Parsec.Token hiding (whiteSpace)-}
{-import Text.ParserCombinators.Parsec.Language (emptyDef)    -}
import Control.Applicative hiding ((<|>),(<?>), many)

eol = try (string "\n")
    <|> try (string "\r")
    <|> try (string "\n\r")
    <?> "End of Line"

whiteSpace = string " "
            <|> string "\t"
            <|> eol
            <?> "Whitespace"

space :: GenParser Char st Char
space = oneOf "         "

spaces = many space

spaces1 = many1 space

notChar x = notFollowedBy x *> anyChar

noWhiteSpace = notChar whiteSpace

noEol = notChar eol
